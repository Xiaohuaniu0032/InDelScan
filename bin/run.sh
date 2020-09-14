python3 InsScan.py -bam /Users/lffu/git_repo/SimVar2Fastq/test/test.sort.Ins.30.bam -n test_sample -outdir $PWD >pysam.read.log

less test_sample.softclip.read.txt|awk 'NR>1 {print $1}'|sort|uniq -c|wc -l
