# RSpec for BioRuby-GFF3-Plugin. Run with something like:
#
#   ruby -I ../bioruby/lib/ ~/.gems/bin/spec spec/gff3_assemble_spec.rb 
#
# Copyright (C) 2010 Pjotr Prins <pjotr.prins@thebird.nl>
#
$: << "../lib"

require 'bio/db/gff/gffdb'

include Bio::GFFbrowser

FASTAFILE="test/data/gff/MhA1_Contig1133.fa"  
GFF3FILE="test/data/gff/MhA1_Contig1133.gff3"

describe GFFdb, "Assemble CDS" do
  before :all do 
    gffdb = Bio::GFFbrowser::GFFdb.new(GFF3FILE, :fasta_filename => FASTAFILE)
    @gff = gffdb.assembler
    @gff.parse
    @contigsequence = @gff.sequencelist["MhA1_Contig1133"]
    @cdslist = {}
    @gff.each_CDS do | id, reclist, component |
      @component = component
      @cdslist[id] = reclist
    end
  end

  it "should have the single contig" do 
    @gff.sequencelist.size.should == 1
    @gff.sequencelist["MhA1_Contig1133"].should_not == nil
    @gff.sequencelist["MhA1_Contig1133"].size.should == 33905
  end

  it "should have CDS 7838:7980" do
    recs = @cdslist['cds:MhA1_Contig1133.frz3.gene4']
    cds0 = recs[0]
    cds0.start.should == 7838
    cds0.end.should == 7980
    cds0.frame.should == 0
    cds0.seqname.should == 'MhA1_Contig1133'
  end
  it "should have CDS 8065:8308" do
    recs = @cdslist['cds:MhA1_Contig1133.frz3.gene4']
    cds1 = recs[1]
    cds1.start.should == 8065
    cds1.end.should == 8308
    cds1.frame.should == 1
    cds1.strand.should == '+'
    cds1.seqname.should == 'MhA1_Contig1133'
  end
  it "should translate CDS 7838:7980 (in frame 0)" do
    recs = @cdslist['cds:MhA1_Contig1133.frz3.gene4']
    cds0 = recs[0]
    cds0.seqname.should == 'MhA1_Contig1133'
    seq = @gff.assemble(@contigsequence,@component.start,[cds0])
    seq.should == "TGCGTCCTTTAACAGATGAAGAAACTGAAAAGTTTTTCAAAAAACTTTCAAATTATATTGGTGACAATATTAAACTTTTATTGGAAAGAGAAGATGGAGAATATGTTTTTCGTTTACATAAAGACAGAGTTTATTATTGCA"
    aaseq = @gff.assembleAA(@contigsequence,@component.start,[cds0])
    aaseq.should == "RPLTDEETEKFFKKLSNYIGDNIKLLLEREDGEYVFRLHKDRVYYCX"
  end
  it "should translate CDS 8065:8308 (in frame 1)" do
    recs = @cdslist['cds:MhA1_Contig1133.frz3.gene4']
    cds1 = recs[1]
    seq = @gff.assemble(@contigsequence,@component.start,[cds1])
    seq.should == "GAAAAATTAATGCGACAAGCAGCATGTATTGGACGTAAACAATTGGGATCTTTTGGAACTTGTTTGGGTAAATTCACAAAAGGAGGGTCTTTCTTTCTTCATATAACATCATTGGATTATTTGGCACCTTATGCTTTAGCAAAAATTTGGTTAAAACCACAAGCTGAACAACAATTTTTATATGGAAATAATATTGTTAAATCTGGTGTTGGAAGAATGAGTGAAGGGATTGAAGAAAAACA"
    aaseq = @gff.assembleAA(@contigsequence,@component.start,[cds1])
    # note it should handle the frame shift and direction!
    aaseq.should == "EKLMRQAACIGRKQLGSFGTCLGKFTKGGSFFLHITSLDYLAPYALAKIWLKPQAEQQFLYGNNIVKSGVGRMSEGIEEKX"
    # we are going to force a change of direction
    cds1rev = cds1
    cds1rev.strand = '-'
    seq = @gff.assemble(@contigsequence,@component.start,[cds1rev])
    seq[0..3].should == "ACAA"
    aaseq = @gff.assembleAA(@contigsequence,@component.start,[cds1rev])
    aaseq.should_not == "EKLMRQAACIGRKQLGSFGTCLGKFTKGGSFFLHITSLDYLAPYALAKIWLKPQAEQQFLYGNNIVKSGVGRMSEGIEEKX"
  end
  it "should assemble 3 CDSs for MhA1_Contig1133.frz3.gene4"
end

