import glob
from pathlib import Path
import os

import argparse
from subprocess import check_output
import eos_posteriors 


executable_file = "/home/weizmann.kiendrebeogo/NMMA/Anna-EOS/eos_posteriors.sh"

def main():

    parser = argparse.ArgumentParser(
        description="Determine EOS_posteriors from Lamda_tilde."
    )
    parser.add_argument(
        "--outdir", 
        type=str, 
        default="outdir",
        help="Path to the output directory"
    )
    parser.add_argument(
        "--GWsamples",
        type=str,
        required=True,
        help="Data from the posterior of the GW injection",
    )
    parser.add_argument(
        "--EOSpath", 
        type=str, 
        required=True, 
        help="The EOS data"
    )
    parser.add_argument(
        "--condor-dag-file",
        type=str,
        required=True,
        help="The condor dag file to be created"
    )
    parser.add_argument(
        "--condor-sub-file",
        type=str,
        required=True,
        help="The condor sub file to be created"
    )
    parser.add_argument(
        "--bash-file", 
        type=str, 
        required=True, 
        help="The bash file to be created"
    )
    args = parser.parse_args()


    outdir = args.outdir            
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
        
    logdir = os.path.join(outdir, f"logs")  
    if not os.path.isdir(logdir):
        os.makedirs(logdir)

    path = Path(args.GWsamples) 
    GWsamples = glob.glob(str(path/'**/*.json'), recursive=True)
    
    number_jobs = len(GWsamples)
    job_number = 0
    
    fid =  open(args.condor_dag_file, "w")
    fid1 = open(args.bash_file, "w")
    
    for gwsamples in GWsamples[0:2]:
                  
        fid.write("JOB %d %s\n" % (job_number, args.condor_sub_file))
        fid.write("RETRY %d 3\n" % (job_number))
        fid.write(
            'VARS %d jobNumber="%d"  outdir="%s" GWsamples="%s"\n'
            % (job_number, job_number, outdir,  gwsamples))
        fid.write("\n\n")
        job_number = job_number + 1
        
        fid1.write(f"{executable_file} --outdir {outdir} --GWsamples {gwsamples} --EOSpath {args.EOSpath}\n".format()
        )
        
    fid.close()
    fid1.close()
    
    fid = open(args.condor_sub_file, "w")
    fid.write("executable = %s\n" % (executable_file))
    fid.write(f"output = {logdir}/out.$(jobNumber)\n")
    fid.write(f"error =  {logdir}/err.$(jobNumber)\n")

    fid.write(
        f"arguments = --outdir $(outdir)  --GWsamples $(GWsamples) --EOSpath {args.EOSpath}\n".format()
        )
    fid.write('requirements = OpSys == "LINUX"\n')
    fid.write("request_memory = 8192\n")
    fid.write("request_disk = 500 MB\n")
    fid.write("request_cpus = 1\n")
    fid.write("accounting_group = ligo.dev.o2.burst.allsky.stamp\n")
    fid.write("notification = nevers\n")
    fid.write("getenv = true\n")
    fid.write("log = /local/%s/eos_to_lambda_posteriors_condor.log\n" % os.environ["USER"])
    fid.write("+MaxHours = 24\n")
    fid.write("universe = vanilla\n")
    fid.write("queue 1\n")


if __name__ == "__main__":
    main()
