#!/bin/sh
echo "subsample sync file"
perl /software/x86_64/popoolation2/subsample-synchronized.pl --input <(bzcat ../30.popoolation.sync/Bolt_Non_java.sync.bz2) --output subsample.sync --target-coverage 50 --max-coverage 1000 --method withreplace
