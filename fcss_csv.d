import std.regex, std.stdio, std.algorithm, std.range;

void main(){
    auto lines = stdin.byLine.map!(line=>line.idup).array;
    foreach(i, line; lines){
        auto m = line.matchFirst(`(\d+.mol) (\d+)`);
        if(m){
            auto codes = lines[i+1].splitter(`\s+`).array;
            write(m[1], ';', m[2]);
            foreach(c; codes) write(';', c);
            writeln();
        }
    }
}