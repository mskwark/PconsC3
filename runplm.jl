#!/home/skwarkmj/sw/julia/julia

using ArgParse
using GaussDCA
using PlmDCA


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--gaps", "-g"
            help = "Fraction of gaps to use"
        "--di"
            help = "an option without argument, i.e. a flag"
            action = :store_true
        "arg1"
            help = "input alignment"
            required = true
        "arg2"
            help = "lambdaJ 0"
            required = true
        "arg3"  
            help = "Output file"
            required = true
    end

    return parse_args(s)
end

function main()
   parsed_args = parse_commandline()
   x = PlmDCA.plmdca(parsed_args["arg1"], lambdaJ=float(parsed_args["arg2"]))
   score=x.score

   io = open(parsed_args["arg3"], "w")
   for s in score
        @printf(io, "%d,%d,%f\n",s[1],s[2],s[3])
   end
   close(io)
end
  

main()
