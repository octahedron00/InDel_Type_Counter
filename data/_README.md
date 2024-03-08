# Data folder of InDel_Type_Counter:

To run the program, you must put the NGS data in here:

only file with name '.fastq' and '.fastq.gz' will be read by the program,

and other files like this will not be used as a data.

### main.py -x R2 -x undetermined 

To ignore some files after getting all data from the NGS raw result,

you can delete some of it manually, or using the config in the main.py code.

the config [-x R2] will ignore files with name 'R2', which mainly means Read 2.
