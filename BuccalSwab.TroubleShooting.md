# How many forwards reads are there?
I've tried to adapt Ludo's concatanation code, to allow me to extract the number of forwards reads for each sample. It was relatively hard, I've never used python.....

Some notes:
- os.popen(...) allows me to run shell commands within the python environment
- f"..." are f-strings which means python with substitue the {...} for actual variables we call in the python script
- with ... as fwd_count this allows us to cnnect between the shell command and the rest of our python script
- : is used to sginal the start of a block of code in python (like in if or else statement), they are not for single instructions
- with is another block structure like if or else, it allows us to use a file and then close it automitcally even when errors occur
  - you have to specifically close files in python (?)
- indentations are used in python to define code blocks, its very confusing

```bash
cd /home/mulha552/uoo04306/frogs_gbs
nano forwardreads.py
```
Copy the code below into 'forwardsreads.py
Note: the correct way to count the number of forwards reads would be to count the total number of lines and divide by four, because @ is a PHred score we are potentially counting additional lines. However, this is unlikely to matter for our purposes. 
```python
#!/usr/bin/env python3

import os
os.chdir('/home/mulha552/uoo04306/frogs_gbs') #python for 'cd'
allfiles_to_count = os.listdir("samples/")
with open("forward_counts.txt", "w") as output_file:
    with open("popmap.txt") as f:
        for line in f:
            samplename=line.split("\t")[0]
            print (samplename) # should be commented out in a slurm job, but I can use this as a check things are running.
            checkfiles=[filename for filename in allfiles_to_count  if filename.startswith(samplename)]
            if len (checkfiles)!=4: # that was weirdly complicated because some sample name are contained in others different ways, but the vcheck above solve it uysing the rem file
                print(line)
                raise Exception
          
            forward_counts = 0
            for checkfile in checkfiles:
                if checkfile.endswith(".1.fq.gz"):
                    file_path = os.path.join("samples", checkfile)
                    if os.path.exists(file_path):
                        with os.popen(f"zcat {file_path} | grep -c '^@'") as fwd_count:
                            forward_counts += int(fwd_count.read().strip())
                    else:
                        print(f"File {file_path} doesn't exist.")
            
            output_file.write(f"{samplename}\t{forward_counts}\n")

```
Now we'll run it
```bash
module load Python #need to load a python environment//ipython console
python3 forwardreads.py
```
# How many forwards reads are there if I demultiplex, without cutadapt?
Ok, now we know there is no obvious difference in the number of forwards reads following our orignal pipeline that would lead to the difference in coverage depth; We now need to understand whether this difference is because of lab or sequencing errors or if it is genuine. To do this I will demultiplex the raw.reads without using cutadapt (or any cleaning option like -c) if there is a difference we know this is a sequencing or lab error (maybe during normalization ?) if not then buccal swabs would have to have a lower percentage of reads lost during CutAdapt ..... I beleive the latter is more likely. 

Run as a SLURM Job. 
```bash
#mkdir samples_noQC
cd /home/mulha552/uoo04306/frogs_gbs
module load Stacks #2.61
process_radtags -P   -p raw/ -o samples_noQC/ -b ../barcodes.txt -e pstI -r  --inline-inline
```
```bash
sh deplex_noQC.sl # make sure it works
squeue -u mulha552 # check job is running
```
Once complete modify "forwardreads.py" to count the number of forwards reads without quality control.

