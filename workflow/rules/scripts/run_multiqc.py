import subprocess
import argparse
import glob
import os
    
def generate_file_list(outdir: str, file_list_filename: str):
    """Generate file list for MultiQC"""
    multiqc_dir = os.path.join(outdir, "multiqc")
    file_list = os.path.join(multiqc_dir, file_list_filename)
    
    with open(file_list, "w") as f:
        analysis_dir = [
            (os.path.join(outdir, "rsem_quant"), "**/*.cnt"),
            (os.path.join(outdir, "align"), "**/*_genome_flagstat.txt"),
            (os.path.join(outdir, "processed"), "*.json"),
            (os.path.join(outdir, "align"), "**/*_Log.final.out")
        ]
        for dir, pattern in analysis_dir:
            files = glob.glob(os.path.join(dir, pattern), recursive=True)
            print("\n".join(files))
            for file in files:
                f.write(file + "\n")
    
def run_multiqc(multiqc_config: str, outdir: str, file_list_filename: str, overwrite: bool):
    """Run MultiQC"""
    multiqc_dir = os.path.join(outdir, "multiqc")
    file_list = os.path.join(multiqc_dir, file_list_filename)
    cmd = ["multiqc", "--config", multiqc_config, "--file-list", file_list, "--outdir", multiqc_dir]
    if overwrite:
        cmd.append("--force")
    subprocess.check_call(cmd)

def main(args):
    multiqc_dir = os.path.join(args.outdir, "multiqc")
    os.makedirs(multiqc_dir, exist_ok=True)
    
    generate_file_list(outdir=args.outdir, file_list_filename=args.file_list)
    run_multiqc(multiqc_config=args.multiqc_config, outdir=args.outdir, file_list_filename=args.file_list, overwrite=args.overwrite)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Prepare custom QC data (*_mqc.yaml) and file list for MultiQC.")
    parser.add_argument("--multiqc-config", help="Path to the MultiQC config file")
    parser.add_argument("--outdir", default="results", help="Path to the output directory")
    parser.add_argument("--file-list", required=False, default="file_list.txt", help="Name of the file list file")
    parser.add_argument("--overwrite", default=False, action="store_true", help="Overwrite any existing multiqc data")
    args = parser.parse_args()
    main(args)