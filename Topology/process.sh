#!/bin/bash
INPUT_FILE="$1"
P="${2:-8}"

fetch_line() {
    line="$1"
    pdb=$(echo "$line" | cut -d ";" -f2)

    b=$(curl -s "https://data.rcsb.org/rest/v1/core/entry/${pdb}")

    resolution=$(echo "$b" | awk '{
        a=index($0, "\"resolution_combined\":"); 
        if(a>0){
            rcl=length("\"resolution_combined\":"); 
            c=substr($0,a+rcl+1); 
            b=index(c,","); 
            print substr($0,a+rcl+1,b-2);
        }
    }')

    method=$(echo "$b" | awk '{
        a=index($0,"\"experimental_method\":\"");
        if(a>0){
            rcl=length("\"experimental_method\":\"");
            c=substr($0,a+rcl);
            b=index(c,"\"");
            print substr(c,1,b-1);
        }
    }')

    echo "${line};${resolution};${method}"
}

export -f fetch_line

cat "$INPUT_FILE" | tail -n +1 | xargs -I {} -P "$P" bash -c 'fetch_line "$@"' _ {}
