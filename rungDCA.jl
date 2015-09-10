#!/home/mjs/sw/julia/julia

using ArgParse
using GaussDCA

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
            help = "output file"
            required = true
    end

    return parse_args(s)
end

function main()
   parsed_args = parse_commandline()
   if parsed_args["di"]
	scores = gDCA(parsed_args["arg1"], min_separation=1, score = :DI)
   else
	scores = gDCA(parsed_args["arg1"], min_separation=1)
   end
   printrank(parsed_args["arg2"], scores)
end

main()
