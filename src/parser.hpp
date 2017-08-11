#pragma once
#include <string>
#include <sstream>
#include <stdexcept>
#include <iostream>
#include <ctype.h>
#include <cmath>

using std::pow;
using std::string;

void error(const string& msg);
void warning(const string& msg);

struct FormatSpec
{
    int n, m;
};

const char* fetchSpec(const char* fmt, FormatSpec& spec);

struct Parser
{
    Parser(std::istream& in) :
    inp(in)
    {
        next();
        lineCnt = 0;
    }

    // read a line and fit prefix to 'fmt'
    void matchfln(const char* fmt)
    {
        while (*fmt){
            if (front() != *fmt)
                error("format error - expected " + string(1, *fmt));
            fmt++;
            next();
        }
        skipLine();
    }

    // read a line and fit prefix to 'fmt' poppulating values of args
    template<class F, class... T>
    void matchfln(const char* fmt, F& a, T&... args)
    {
        while (*fmt)
        {
            FormatSpec spec;
            fmt = fetchSpec(fmt, spec);
            if (spec.n) //fetched spec
            {
                if (spec.m)
                    parse(a, spec.n, spec.m);
                else
                    parse(a, spec.n);
                return matchfln(fmt, args...);
            }
            char marker = *fmt;
            if (front() != marker)
                error("format error - pattern mismatch expected:" + string(1, front()));
            fmt++;
            next();
        }
    }

    string quotedString(char delim = ',')
    {
        string field;
        if (front() == '"')
        {
            next();
            field = readUpTo('"');
            if (front() != delim)
                error(string("expected '") + delim + "' not '"+front()+"'");
            next();
        }
        else
            field = readUpTo(delim);
        return field;
    }

    string line()
    {
        return readUpTo('\n');
    }

    bool eof()
    {
        return front_ == -1;
    }

private:
    std::istream& inp;
    int front_;
    int lineCnt;

    char front()
    {
        if (front_ == -1)
            errorEof();
        return front_;
    }

    void error(const string& msg)
    {
        std::stringstream s;
        s << "Line " << lineCnt << ":" << msg;
        throw std::logic_error(s.str());
    }

    void errorEof()
    {
        error("unexpected end of file");
    }

    void next()
    {
        front_ = inp.get();
    }

    int skipSpace(int maxSpace)
    {
        for (int i = 0; i < maxSpace; i++)
        {
            if (!isspace(front()))
                return i;
            next();
        }
        return maxSpace;
    }

    void skipLine()
    {
        while (front() != '\n')
            next();
        lineCnt++;
        next();
    }

    string readUpTo(char delim)
    {
        string line;
        while (!eof() && front() != delim)
        {
            line.push_back(front());
            next();
        }
        if (delim == '\n')
        {
            lineCnt++;
            //strip \r if any
            if(line.size() > 1 && line[line.size()-1] == '\r')
                line.resize(line.size()-1);
        }
        if (!eof())
            next();
        return line;
    }

    template<class T>
    void parse(T&, int)
    {
        error("unsupported type combination");
    }

    template<class T>
    void parse(T&, int, int)
    {
        error("unsupported type combination");
    }

    void parse(string& s, int n)
    {
        s.clear();
        int k = skipSpace(n);
        for (int i = k; i < n; i++)
        {
            s.push_back(front());
            next();
        }
    }

    template<int N>
    void parse(char(&arr)[N], int n)
    {
        if (n + 1 > N)
            error("format - wrongly sized array");
        int k = skipSpace(n);
        for (int i = k; i < n; i++)
        {
            arr[i - k] = front();
            next();
        }
        arr[n - k] = 0;
    }

    void parse(int& target, int n)
    {
        int value = 0;
        for (int i = 0; i < n; i++)
        {
            if (isdigit(front()))
            {
                value = front() - '0' + value * 10;
            }
            else if (!isspace(front())) //may skip space as well            
                error("format error - expected digit:" + string(1, front()));
            else if (front() == '\n')
                break; //HACK for bad formats
            next();
        }
        target = value;
    }

    void parse(double& target, int n, int m)
    {
        char buf[64];
        int k = n + m + 1;
        //assert(n + m + 2 < 64);
        for (int i = 0; i < k; i++)
        {
            buf[i] = front();
            next();
        }
        buf[k] = 0;
        target = atof(buf);
        //target = negative ? -ret : ret;
    }

};
