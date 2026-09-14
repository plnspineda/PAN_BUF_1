echo -e "Sample\tRealTime_sec\tCPU_sec"

find work -name ".command.log" | while read log; do

    main=$(grep "Real time:" "$log")
    [ -z "$main" ] && continue

    dir=$(dirname "$log")
    sample=$(echo "$main" | grep -oE '(ERR|SRR)[0-9]+' | head -1)

    realtime=$(echo "$main" | sed -E 's/.*Real time: ([0-9.]+) sec;.*/\1/')
    cpu=$(echo "$main" | sed -E 's/.*CPU: ([0-9.]+) sec/\1/')

    echo -e "$sample\t$realtime\t$cpu"

done
