import subprocess
import argparse
import vcf
import os
import shutil
import logging 

__author__ = "Vinay Kaikala"

class GenerateSymlink():
    "Generates symlinks required for IGVnav inputs"
    
    def __init__(self, outputdirname):
        """ 
        :params: outputdirname -  sample output directory name 
        """
        self.outputdirname = outputdirname

    def generateIGVsymlink(self, *args):
        """
        Generates the symlinks for give list of tuples
        :params: args - ('directoryname', 'patters')
        """   
        igvnav_dirname_dst = os.path.join(self.outputdirname, 'IGVnav')
        src_dir = os.path.abspath(self.outputdirname)
        logging.info("Generating IGVNav Symlinks in : " + igvnav_dirname_dst ) 
        if not os.path.exists(igvnav_dirname_dst): os.mkdir(igvnav_dirname_dst)
        try:
            symlinks = (('variants','.vep.vcf'),('bams','-nodups.bam'), ('bams','.overlapped.bam'), ('variants','.vep.vcf'), \
                      ('svs/igv','.mut'), ('svs','.gtf'),('svs','.bam'),('svs/svaba', '.contigs.bam')) +  args
            for each_input in symlinks:
                dir_name = os.path.join(src_dir,each_input[0])
                self.create_symlink(dir_name, src_dir, igvnav_dirname_dst, each_input[1])
            logging.info("Created IGVNav Symlinks : " + igvnav_dirname_dst ) 
        except Exception as e:
            logging.info(e)
        return  

    def create_symlink(self, travers_dir_name, src_dir, igvnav_dirname_dst, suffix):
        "Recursively Traverse through the directory and create symlink"
        for root, dirs, files in os.walk(travers_dir_name):
            for each_file in files:
                if each_file.endswith(suffix) and not os.path.exists(os.path.join(igvnav_dirname_dst,each_file)):
                    try:		
                        os.symlink(os.path.join(root,each_file), os.path.join(igvnav_dirname_dst,each_file))
                    except Exception as e:
                        logging.info(e) 
        return

