# This has been tested on Ubuntu 16.04 and CentOS 6 and works optimally on Ubuntu, and does the trick on CentOS, but with some minor visual differences.
HASH=$1
for D in `find cromwell-executions/*/$HASH/ -mindepth 1 -maxdepth 1 -type d -printf "%Tx %.9TX %p\n" | sort | awk {'print $6'}`; do \
START=$(find $D -maxdepth 2 \( -name 'shard*' -o -name 'call*' \) -printf '%Tx at %.10TX\n' | sort | head -n 1) && \
echo "$(basename $D | cut -c 6-) started on $START" && \
STOP=$(find $D -maxdepth 3 -type f \( -name "*.sam" -name "*.bam" -o -name "*.vcf" -o -name "stderr" -o -name "stdout" \) -printf 'And stopped on %Tx at %.9TX \n' | sort | tail -n 1) && \
echo $STOP && \
SIZE=$(find $D -maxdepth 4 -type f \( -name "*.recal" -o -name "*.tranches" -o -name '*.*am' -o -name '*.vcf' -o -name '*.grp' \) -exec du -sh {} \; | sort -h | awk {'print $1'}) && \
echo $SIZE && \
SHARDS=$(find $D -mindepth 1 -maxdepth 2 -type d -name "shard*" | wc -l) && \

	if [ $SHARDS -gt 0 ]; then \
		echo "$SHARDS shards created"; \
	elif [ $SHARDS -eq  0 ]; then \
		:
	fi && \

echo ""; done 
# et voila, run this with "watch -d sh scriptname.sh" to get automatic updates every 2 seconds.
# It would be nice to figure out a way to calculate execution time per tool as well...

# I couldn't put the comments on the same line because it broke the script, I'll have to improve this somehow
# for loop that finds all "call-toolname" folders
# the creation time for each "call-toolname" directory is stored here
# the cut command deletes the "call-" part and the rest of the command prints a nice message in the form of "toolname started on Y/M/D at H:M:S"
# the newest file creation time is used to determine when the tool last updated or finished its execution
# this just prints a message saying when the tool finished/stopped/updated a file in the form of "And stopped on Y/M/D at H:M:S"
# this line counts the number of created shards in scatter/gather processes, and the if/else statement below makes sure to not print "0 shards created" for tools that don't have any shards
