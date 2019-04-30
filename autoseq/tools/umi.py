import logging
import sys
import uuid

from pypedream.job import Job, repeat, required, optional, conditional

class FastqToBam(Job):
	def __init__(self):
		Job.__init__(self)
		self.input_fastq1 = None
		self.input_fastq2 = None
		self.output_bam = None
		self.sample = None
		self.library = None

	def command(self):

		tmpdir = "{}/{}".format(self.scratch, uuid.uuid4())

		cmd = "fgbio -Xmx10g -XX:+AggressiveOpts -XX:+AggressiveHeap -XX:ParallelGCThreads=8 --tmp-dir {} ".format(tmpdir) + \
			" FastqToBam " + \
			" -i {}  {} ".format(self.input_fastq1, self.input_fastq2)  + \
			" -o " + self.output_bam + \
			" --sample {} --library {} ".format(self.sample, self.library) + \
			" -r 3M2S+T 3M2S+T -s true "

		rm_tmpdir = "rm -rf {} ".format(tmpdir)

		return cmd + " && " + rm_tmpdir 

class AlignUnmappedBam(Job):
	def __init__(self):
		Job.__init__(self)
		self.input_bam = None
		self.reference_genome = None
		self.threads = None
		self.output_bam = None

	def command(self):

		tmpdir = "{}/{}".format(self.scratch, uuid.uuid4())

		picard_cmd = "picard SamToFastq I={} F=/dev/stdout INTERLEAVE=true TMP_DIR={} ".format(self.input_bam, tmpdir) 

		bwa_cmd = "bwa mem -p -t {} {} /dev/stdin ".format(self.threads, self.reference_genome)

		picard_merge_cmd = "picard -Xmx10g MergeBamAlignment UNMAPPED={} ALIGNED=/dev/stdin".format(self.input_bam) + \
							" O={} ".format(self.output_bam) + \
							" R={} ".format(self.reference_genome) + \
							" SO=coordinate ALIGNER_PROPER_PAIR_FLAGS=true MAX_GAPS=-1 ORIENTATIONS=FR CREATE_INDEX=true " + \
							" TMP_DIR=".format(tmpdir)

		rm_tmpdir = "rm -rf {} ".format(tmpdir)

		return " | ".join([picard_cmd, bwa_cmd, picard_merge_cmd]) + " && " + rm_tmpdir 

class GroupReadsByUmi(Job):
	def __init__(self):
		Job.__init__(self)
		self.input_bam = None
		self.output_bam = None
		self.output_histogram = None
		

	def command(self):

		tmpdir = "{}/{}".format(self.scratch, uuid.uuid4())

		cmd = "fgbio -Xmx10g -XX:+AggressiveOpts -XX:+AggressiveHeap -XX:ParallelGCThreads=8 --tmp-dir {} GroupReadsByUmi ".format(tmpdir) + \
			  " -i " + self.input_bam + \
			  " -o " + self.output_bam + \
			  " --strategy paired --family-size-histogram " + self.output_histogram

		rm_tmpdir = "rm -rf {} ".format(tmpdir)

		return cmd + " && " + rm_tmpdir 

class CallDuplexConsensusReads(Job):
	def __init__(self):
		Job.__init__(self)
		self.input_bam = None
		self.output_bam = None

	def command(self):
		tmpdir = "{}/{}".format(self.scratch, uuid.uuid4())

		cmd = "fgbio -Xmx10g -XX:+AggressiveOpts -XX:+AggressiveHeap -XX:ParallelGCThreads=8 --tmp-dir {} CallDuplexConsensusReads".format(tmpdir) + \
			  " -i " + self.input_bam + \
			  " -o " + self.output_bam + \
			  " --min-reads 1 1 0 --min-input-base-quality 30 "

		rm_tmpdir = "rm -rf {} ".format(tmpdir)

		return cmd + " && " + rm_tmpdir 

class FilterConsensusReads(Job):
	def __init__(self):
		Job.__init__(self)
		self.input_bam = None
		self.reference_genome = None
		self.output_bam = None

	def command(self):
		tmpdir = "{}/{}".format(self.scratch, uuid.uuid4())

		cmd = "fgbio -Xmx10g -XX:+AggressiveOpts -XX:+AggressiveHeap -XX:ParallelGCThreads=8 --tmp-dir {} FilterConsensusReads".format(tmpdir) + \
		      " -i {}".format(self.input_bam) + \
		      " -o " + self.output_bam + \
		      " --ref " + self.reference_genome + \
		      " --min-reads 1 1 0 --reverse-per-base-tags true --max-read-error-rate 1 " + \
		      " --max-base-error-rate 0 --min-base-quality 30 --require-single-strand-agreement true" 

		rm_tmpdir = "rm -rf {} ".format(tmpdir)

		return cmd + " && " + rm_tmpdir 
   
class ClipBam(Job):
	def __init__(self):
		Job.__init__(self)
		self.input_bam = None
		self.output_bam = None
		self.metrics_txt = None
		self.reference_genome = None
	
	def command(self):
		tmpdir = "{}/{}".format(self.scratch, uuid.uuid4())

		cmd = "fgbio -Xmx10g -XX:+AggressiveOpts -XX:+AggressiveHeap -XX:ParallelGCThreads=8 --tmp-dir {} ClipBam ".format(tmpdir) + \
			  " -i " + self.input_bam + \
			  " -o " + self.output_bam + \
			  " -m " + self.metrics_txt + \
			  " --ref " + self.reference_genome + \
			  " --clip-overlapping-reads true "

		rm_tmpdir = "rm -rf {} ".format(tmpdir)

		return cmd + " && " + rm_tmpdir  
