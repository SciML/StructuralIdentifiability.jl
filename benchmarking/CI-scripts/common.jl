function dump_results(file, key)
    return open(file, "w") do out
        println(out, key)
        for problem in suite
            problem_name = problem.problem_name
            result = problem.result
            type = problem.type
            println(out, "$problem_name:$(String(type)):$(join(map(string, result), ","))")
        end
    end
end
