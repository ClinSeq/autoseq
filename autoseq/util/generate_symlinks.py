import subprocess
import argparse
#import vcf
import os
import shutil
import logging
import re
from collections import defaultdict

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
            symlinks = (('variants','.vep.vcf'),('bams','-nodups.bam'), ('bams','-nodups.bam.bai'),
                        ('bams','.overlapped.bam'), ('bams','.overlapped.bai'), ('variants','.vep.vcf'), ('cnv', '.bedGraph'), ('variants', '.bedGraph'),
                      ('svs/igv','.mut'), ('svs','.gtf'), ('svs','.bam'), ('svs','.bai'), ('', 'igvnav-input.txt')) + args
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

    def get_all_files(self, dir_path):
        igv_session_file_list = defaultdict(lambda : defaultdict(list))
        all_files = [('bam_common', 'bam_nodups', '.*-nodups.bam$'),
                     ('sv', 'bam_cfdna', '.*-CFDNA-.*(DEL|DUP|INV|TRA|contigs.sort).bam$'),
                     ('sv', 'bam_normal', '^(?:(?!CFDNA).)*(DEL|DUP|INV|TRA|contigs.sort).bam$'),
                     ('sv', 'mut_svict_cfdna', '.*-CFDNA-.*.svict_SR8.mut$'),
                     ('sv', 'mut_svict_normal', '^(?:(?!CFDNA).)*.svict_SR8.mut$'),
                     ('sv', 'mut_sava','.*_(germline|somatic)_svaba.mut$'),
                     ('sv', 'mut_svcaller_cfdna', '.*-CFDNA-.*_svcaller.mut$'),
                     ('sv', 'mut_svcaller_normal', '^(?:(?!CFDNA).)*_svcaller.mut$'),
                     ('sv', 'mut_lumpy', '.*lumpy_len500_SU24.mut$'),
                     ('sv', 'gtf_cfdna', '.*-CFDNA-.*(DEL|DUP|INV|TRA).gtf$'),
                     ('sv', 'gtf_normal', '^(?:(?!CFDNA).)*(DEL|DUP|INV|TRA).gtf$'),
                     ('snps', 'bam_cfdna', '.*-CFDNA-.*.clip.overlapped.bam$'),
                     ('snps', 'bam_normal', '^(?:(?!CFDNA).)*clip.overlapped.bam$'),
                     ('snps', 'vep', '.*.all.(somatic|germline).vep.vcf$'),
                     ('cnv', 'flank_profile_cfdna', '.*-CFDNA-.*_profile.bedGraph'),
                     ('cnv', 'flank_profile_normal', '^(?:(?!CFDNA).)*_profile.bedGraph'),
                     ('cnv', 'flank_cnv_cfdna', '^(?:(?!CFDNA).)*_segments.bedGraph'),
                     ('cnv', 'flank_cnv_normal', '.*-CFDNA-.*_segments.bedGraph'),
                     ('asf', 'flank_asf', '.*germline-variants-somatic-afs.bedGraph'),

                     ]

        for variant_type, file_type, regex in all_files:
            files_list = list(filter( lambda x: re.match(regex, x) and not x.startswith('.') and not x.endswith('.out'), os.listdir(dir_path)))
            files_list.sort()
            igv_session_file_list[variant_type][file_type]=files_list

        return igv_session_file_list

    def create_igv_session_file(self, igv_session_master="/nfs/PROBIO/for_igv/igv_session_master.xml", igv_session_sv_master="/nfs/PROBIO/for_igv/igv_session_sv_master.xml"):
        """
        Create IGV seesion file in xml format for given sample
        """
        resource_path = """
        <Resource path="{path_name}"/>
        """
        panel_str = """
        <Panel height="{panel_height}"  width="{panel_width}" name="{panel_name}">
            {tracks}
        </Panel>
        """
        tracks_bam_str= """
        <Track autoScale="true" clazz="org.broad.igv.sam.CoverageTrack" color="175,175,175" colorScale="ContinuousColorScale;0.0;927.0;255,255,255;175,175,175" fontSize="10" id="{bam_file_full}_coverage" name="{bam_file} Coverage" snpThreshold="0.0001" visible="true">
            <DataRange baseline="0.0" drawBaseline="true" flipAxis="false" maximum="1503.0" minimum="0.0" type="LINEAR"/>
        </Track>
        <Track clazz="org.broad.igv.sam.SpliceJunctionTrack" fontSize="10" height="60" id="{bam_file_full}_junctions" name="{bam_file} Junctions" visible="false"/>
        <Track clazz="org.broad.igv.sam.AlignmentTrack" displayMode="EXPANDED" experimentType="OTHER" fontSize="10" id="{bam_file_full}" name="{bam_file}" visible="true">
            <RenderOptions/>
        </Track>  
        """
        sv_gtf_track_str="""
            <Track clazz="org.broad.igv.track.FeatureTrack" color="0,0,178" fontSize="10" id="{gtf_file_full_path}" name="{gtf_file_path}" visible="true"/>
        """
        sv_mut_track_str="""
                <Track clazz="org.broad.igv.track.MutationTrack" color="0,0,178" colorScale="ContinuousColorScale;0.0;12.0;255,255,255;0,0,178" fontSize="10" height="15" id="{mut_file_full_path}_{id_field}" name="{id_field}" visible="true"/>
        """
        snp_vcf = """
            <Track clazz="org.broad.igv.variant.VariantTrack" color="0,0,178" displayMode="EXPANDED" fontSize="10" id="{vcf_full_file_path}" name="{vcf_file_path}" siteColorMode="ALLELE_FREQUENCY" squishedHeight="1" visible="true"/>
        """
        cnv_flank = """
            <Track autoScale="false" clazz="org.broad.igv.track.DataSourceTrack" displayMode="SQUISHED" fontSize="10" height="80" id="{full_path}" name="CNV {full_path}" color="255,51,51" renderer="BAR_CHART" visible="true" windowFunction="mean">
            <DataRange baseline="0.0" drawBaseline="true" flipAxis="false" maximum="1.0" minimum="-1" type="LINEAR"/>
        </Track>
        """
        cnv_flank_profile = """
            <Track  autoScale="false" clazz="org.broad.igv.track.DataSourceTrack" displayMode="SQUISHED" fontSize="10" height="80" id="{full_path}" name="CNV Profile {full_path}" color="255,51,51" renderer="SCATTER_PLOT" visible="true" windowFunction="mean">
            <DataRange baseline="0.0" drawBaseline="true" flipAxis="false" maximum="1.0" minimum="-1" type="LINEAR"/>
        </Track>

        """
        snp_asf = """
            <Track  clazz="org.broad.igv.track.DataSourceTrack" displayMode="SQUISHED" fontSize="10" height="80" id="{full_path}" name="SNP ASF {full_path}" color="255,51,51" renderer="SCATTER_PLOT" color="255,51,51" visible="true" windowFunction="mean">
            <DataRange baseline="0.5" drawBaseline="true" flipAxis="false" maximum="1.0" minimum="-1" type="LINEAR"/>
        </Track>

        """

        
        try:
            logging.info(" Generating IGV Session File ")
            igvnav_dir = os.path.join(self.outputdirname, 'IGVnav')
            if not os.path.exists(igvnav_dir): raise Exception('IGVnav folder not found..')
            if not os.path.exists(igv_session_master): raise Exception('IGV session master file not found : /nfs/PROBIO/for_igv/igv_session_master.xml ' )
            if not os.path.exists(igv_session_master): raise Exception('IGV session master file not found : /nfs/PROBIO/for_igv/igv_session_sv_master.xml ' )
            igv_session_master_str=open(igv_session_master, 'r').read()
            igv_session_sv_master_str=open(igv_session_sv_master, 'r').read()
            igv_session_files = self.get_all_files(igvnav_dir)

            #session file for snps
            snp_resource = ""
            snp_bam_panel = ""
            snp_vep = ""
            for track_type, each_track in [('snps', 'bam_cfdna'), ('snps', 'bam_normal'),
                     ('cnv', 'flank_cnv_normal'),
                     ('cnv', 'flank_profile_normal'),
                     ('cnv', 'flank_cnv_cfdna'),
                     ('cnv', 'flank_profile_cfdna'),
                     ('asf', 'flank_asf'),
                     ('snps', 'vep')]:
                
                for each_file in igv_session_files[track_type][each_track]:
                    full_path = igvnav_dir + '/' + each_file
                    if not os.path.exists(full_path) or not os.path.getsize(full_path):
                        continue
                    snp_resource += resource_path.format(path_name=full_path)
                    if each_track.startswith('bam'):
                        snp_bam_track = tracks_bam_str.format(bam_file_full=full_path, bam_file=each_file)
                        snp_bam_panel += panel_str.format(panel_height=250, panel_width=1800, panel_name=each_file , tracks=snp_bam_track)
                    if each_track.startswith('flank_cnv'):
                        snp_vep += cnv_flank.format(full_path=full_path)
                    if each_track.startswith('flank_profile'):
                        snp_vep += cnv_flank_profile.format(full_path=full_path)
                    if each_track.startswith('flank_asf'):
                        snp_vep += snp_asf.format(full_path=full_path)
                    if  each_track.startswith('vep'):
                        snp_vep += snp_vcf.format(vcf_full_file_path=full_path, vcf_file_path=each_file)

             
            #session file for structural variants
            all_sv_files=[('bam_common', 'bam_nodups'),('sv', 'bam_cfdna'), ('sv', 'bam_normal'),
                     ('cnv', 'flank_cnv_normal'),
                     ('cnv', 'flank_profile_normal'),
                     ('cnv', 'flank_cnv_cfdna'),
                     ('cnv', 'flank_profile_cfdna'),
                     ('asf', 'flank_asf'),
                     ('sv', 'mut_svict_cfdna'), ('sv', 'mut_svict_normal'),
                     ('sv', 'mut_sava'),  ('sv', 'mut_svcaller_cfdna'),
                     ('sv', 'mut_svcaller_normal'),
                     ('sv', 'mut_lumpy'), ('sv', 'gtf_cfdna'),
                     ('sv', 'gtf_normal')]

            sv_resource= ""
            sv_bam_panel= ""
            sv_mut_track= ""
            sv_gtf_track= ""
            for track_type, each_track in all_sv_files:
                for each_file in igv_session_files[track_type][each_track]:

                    full_path = igvnav_dir + '/' + each_file
                    if not os.path.exists(full_path) :
                        continue

                    if not os.path.getsize(full_path):
                        continue

                    #sv_resource += resource_path.format(path_name=full_path)

                    if each_track.startswith('bam'):
                        sv_resource += resource_path.format(path_name=full_path)
                        sv_bam_track = tracks_bam_str.format(bam_file_full=full_path, bam_file=each_file)
                        sv_bam_panel += panel_str.format(panel_height=150, panel_width=2543, panel_name=each_file, tracks=sv_bam_track)

                    if each_track.startswith('mut'):
                        # get content of SDID column to append to id and name fields in the track
                        f = open(full_path, 'r')
                        line = f.readline() # skip the header line
                        line = f.readline()  # get the first variant line
                        f.close()
                        if line:
                            id_field = line.strip().split()[3]  # use the 4th column (SDID column) if there is a variant
                        else:
                            continue  # skip the file if there are no variants
                        
                        sv_resource += resource_path.format(path_name=full_path)
                        sv_mut_track += sv_mut_track_str.format(mut_file_full_path=full_path, id_field=id_field)

                    if each_track.startswith('flank_cnv'):
                        sv_resource += resource_path.format(path_name=full_path)
                        sv_mut_track += cnv_flank.format(full_path=full_path)
                    if each_track.startswith('flank_profile'):
                        sv_resource += resource_path.format(path_name=full_path)
                        sv_mut_track += cnv_flank_profile.format(full_path=full_path)
                    if each_track.startswith('flank_asf'):
                        sv_resource += resource_path.format(path_name=full_path)
                        sv_mut_track += snp_asf.format(full_path=full_path)

                    if each_track.startswith('gtf'):
                        sv_resource += resource_path.format(path_name=full_path)
                        sv_gtf_track += sv_gtf_track_str.format(gtf_file_full_path=full_path, gtf_file_path=each_file)


            sv_session_data=igv_session_sv_master_str.format(add_resource=sv_resource, add_sv_mut_track=sv_mut_track, add_panel=sv_bam_panel, add_sv_gtf_track=sv_gtf_track)
            snp_session_data=igv_session_master_str.format(add_resource=snp_resource, add_panel=snp_bam_panel, add_vcf_track=snp_vep)
            with open(igvnav_dir+'/igv_session_snps.xml', 'w') as fw:
                fw.write(snp_session_data)

            with open(igvnav_dir+'/igv_session_sv.xml', 'w') as fw:
                fw.write(sv_session_data)

            logging.info("Created IGV Session Files..")

        except Exception as e:
            print('error', e)
            logging.info("Error While creating session file : " + str(e))






