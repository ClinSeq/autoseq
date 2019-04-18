from pypedream.job import Job, required, optional, conditional


class Svcaller(Job):
    def __init__(self):
        Job.__init__(self)
        self.input_bam = None
        self.event_type = None
        self.output_bam = None
        self.output_gtf = None
        self.reference_sequence = None
        self.scratch = None
        self.jobname = "svcaller-run-all"

    def command(self):
        
        activate_env_cmd = "source activate svcallerenv "

        run_all_cmd = ("svcaller run-all --tmp-dir {scratch} " +
                      "--event-type {event_type} " +
                      "--fasta-filename {reference_seq} " +
                      "--filter-event-overlap --events-gtf {output_gtf} "
                      "--events-bam {output_bam} {input_bam}").format(
                          scratch=self.scratch,
                          event_type=self.event_type,
                          reference_seq=self.reference_sequence,
                          output_gtf=self.output_gtf,
                          output_bam=self.output_bam,
                          input_bam=self.input_bam,
                      )

        deactivate_env_cmd = "source deactivate"

        return "{} && {} && {}".format(
            activate_env_cmd,
            run_all_cmd,
            deactivate_env_cmd,
        )

class Sveffect(Job):
    def __init__(self):
        Job.__init__(self)
        self.input_del_gtf = None
        self.input_dup_gtf = None
        self.input_inv_gtf = None
        self.input_tra_gtf = None
        self.ts_regions = None
        self.ar_regions = None
        self.fusion_regions = None
        self.output_combined_bed = None
        self.output_effects_json = None
        self.jobname = "sveffect"

    def command(self):
        
        activate_env_cmd = "source activate svcallerenv "
        make_bed_cmd = ("sveffect make-bed " +
                       "--del-gtf {del_gtf} " +
                       "--dup-gtf {dup_gtf} " +
                       "--inv-gtf {inv_gtf} " +
                       "--tra-gtf {tra_gtf} " +
                       "{output_combined_bed}").format(
                           del_gtf=self.input_del_gtf,
                           dup_gtf=self.input_dup_gtf,
                           inv_gtf=self.input_inv_gtf,
                           tra_gtf=self.input_tra_gtf,
                           output_combined_bed=self.output_combined_bed
                       )

        predict_cmd = ("sveffect predict " +
                      "--ts-regions {ts_regions} " +
                      "--ar-regions {ar_regions} " +
                      "--fusion-regions {fusion_regions} " +
                      "--effects-filename {output_effects_json} {combined_effects_bed}").format(
                          ts_regions=self.ts_regions,
                          ar_regions=self.ar_regions,
                          fusion_regions=self.fusion_regions,
                          output_effects_json=self.output_effects_json,
                          combined_effects_bed=self.output_combined_bed,
                      )

        deactivate_env_cmd = "source deactivate"

        return "{} && {} && {} && {}".format(
            activate_env_cmd,
            make_bed_cmd,
            predict_cmd,
            deactivate_env_cmd,
        )

class MantaSomaticSV(Job):
    def __init__(self):
        Job.__init__(self)
        self.input_tumor = None
        self.input_normal = None
        self.tumorid = None
        self.normalid = None
        self.reference_sequence = None
        self.target_bed = None
        self.output_dir = None
        self.jobname = "manta-somatic-sv"
        
    def command(self):
        required("", self.input_tumor)
        required("", self.input_normal)
        required("", self.reference_sequence)

        # configuration - exome param is added as default, Need to be removed if experiment is WGS.
        configure_mantasv = ("configManta.py " + 
                                    "--generateEvidenceBam " + 
                                    "--outputContig " + 
                                    "--exome " + 
                                    "--callRegions {target_bed} " + 
                                    "--normalBam {input_normal} " + 
                                    "--tumorBam {input_tumor} " + 
                                    "--referenceFasta {reference_sequence} " + 
                                    "--runDir {output_dir} ").format(
                                          target_bed=self.target_bed, 
                                          input_normal=self.input_normal, 
                                          input_tumor=self.input_tumor, 
                                          reference_sequence=self.reference_sequence, 
                                          output_dir=self.output_dir
                                      ) 

        cmd = configure_mantasv + " && " + self.output_dir+"/runWorkflow.py -m local -j 20"
        return cmd

class SViCT(Job):
  def __init__(self):
    Job.__init__(self)
    self.input_bam = None
    self.reference_sequence = None
    self.output = None
    self.jobname = "svict-sv-calling"
        
  def command(self):

    cmd = ("svict -i {input_bam} -r {reference_sequence} -o {output} ").format(input_bam=self.input_bam, reference_sequence=self.reference_sequence, output=self.output)
    return cmd

class Svaba(Job):
  def __init__(self):
    Job.__init__(self)
    self.input_normal = None
    self.input_tumor = None
    self.target_bed = None
    self.output_sample = None
    self.reference_sequence = None
    self.threads =  None
    self.jobname = "svaba-sv-calling"

  def command(self):

    cmd = ("svaba run -t {tumor} -n {normal} -G {reference_sequence} -p {threads} -a {output_sample} ").format(
              tumor = self.input_tumor,
              normal = self.input_normal,
              reference_sequence = self.reference_sequence,
              threads = self.threads,
              output_sample = self.output_sample
              ) 

    return cmd
    
class Lumpy(Job):
  def __init__(self):
    Job.__init__(self)
    self.input_normal = None
    self.input_tumor = None
    self.normal_discordants = None
    self.tumor_discordants = None
    self.normal_splitters = None
    self.tumor_splitters =  None
    self.output = None
    self.threads = None
    self.jobname = "lumpy-sv-calling"

  def command(self):

    discordant_cmd = ("samtools view -@ {threads} -b -F 1294 {n_bam} > {n_discordants} " + \
                    " && samtools view -@ {threads} -b -F 1294 {t_bam} > {t_discordants}").format(
                          threads = self.threads,
                          n_bam = self.input_normal,
                          t_bam = self.input_tumor,
                          n_discordants = self.normal_discordants,
                          t_discordants = self.tumor_discordants)

    splitter_cmd = ("samtools view -@ {threads} -h {n_bam} " + \
                   "| extractSplitReads_BwaMem -i stdin | samtools view -@ {threads} -Sb - " + \
                   " > {n_splitters} && samtools view -@ {threads} -h {t_bam}" + \
                   "| extractSplitReads_BwaMem -i stdin | samtools view -@ {threads} -Sb - " + \
                   " > {t_splitters} ").format(
                          threads = self.threads,
                          n_bam = self.input_normal,
                          t_bam = self.input_tumor,
                          n_splitters = self.normal_splitters,
                          t_splitters = self.tumor_splitters
                   )
    lumpy_cmd = ("lumpyexpress -B {t_bam},{n_bam} -S {t_splitters},{n_splitters} " + \
                " -D {t_discordants},{n_discordants} -o {output}").format(
                          threads = self.threads,
                          n_bam = self.input_normal,
                          t_bam = self.input_tumor,
                          n_discordants = self.normal_discordants,
                          t_discordants = self.tumor_discordants,
                          n_splitters = self.normal_splitters,
                          t_splitters = self.tumor_splitters,
                          output = self.output)

    return " && ".join([discordant_cmd, splitter_cmd, lumpy_cmd])
    


    