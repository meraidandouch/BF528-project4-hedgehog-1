#!/bin/bash
for acc in SRR3879604 SRR3879605 SRR3879606; do
    prefetch --max-size 50G $acc &
done
wait  # Wait for all background jobs to complete
