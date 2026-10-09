#! /usr/bin/awk -f
# Split query-sorted SAM records without splitting a read-name group.
BEGIN { i = 1; count = 0 }
{
    if (count >= fragment_size && $1 != previous) {
        i++
        count = 0
    }
    print >> directory "/" i "_" sample ".sam"
    previous = $1
    count++
}
