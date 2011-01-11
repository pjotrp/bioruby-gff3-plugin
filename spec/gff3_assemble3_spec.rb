# RSpec for BioRuby-GFF3-Plugin. Run with something like:
#
#   rspec -I ../bioruby/lib/ spec/gff3_assemble3_spec.rb 
#
# Copyright (C) 2010,2011 Pjotr Prins <pjotr.prins@thebird.nl>
#
$: << "../lib"

require 'bio-gff3'

include Bio::GFFbrowser

GFF3FILE3="test/data/gff/test-cds.gff3"

describe GFF3, "Assemble CDS (extra checks)" do
  before :all do 
    gff3 = Bio::GFFbrowser::GFF3.new(GFF3FILE3)
    @gff = gff3.assembler
    @gff.parse
  end

  it "should translate gene MhA1_Contig1040.frz3.gene29" do
    @contigsequence = @gff.sequencelist["MhA1_Contig1040"]
    @componentlist = {}
    @cdslist = {}
    @gff.each_CDS do | id, reclist, component |
      @componentlist[id] = component
      @cdslist[id] = reclist
    end
    name = "cds:MhA1_Contig1040.frz3.gene"
    recs = @cdslist[name]
    component = @componentlist[name]
    ntseq = @gff.assemble(@contigsequence,component.start,recs,:raw=>true)
    ntseq.should == "TTAATTAATTTGCCTAGAAAAACAAAGGCATAACATGCTTGCAGTCATCATACGGTAAGAGAGAAACCAACGATATGTTAATAATGTTGATGGGGGAATATCCTCATTAGAATTCTTTTTTGGGTGAATTGAAATTGCCATATTATTAGTATTATTAGAAAATATTAAATTTGTTGATAA"
    ntseq = @gff.assemble(@contigsequence,component.start,recs,:codonize=>true)
    ntseq.should == "TTATCAACAAATTTAATATTTTCTAATAATACTAATAATATGGCAATTTCAATTCACCCAAAAAAGAATTCTAATGAGGATATTCCCCCATCAACATTATTAACATATCGTTGGTTTCTCTCTTACCGTATGATGACTGCAAGCATGTTATGCCTTTGTTTTTCTAGGCAAATTAATTAA"
    aaseq = @gff.assembleAA(@contigsequence,component.start,recs)
    aaseq.should == "LSTNLIFSNNTNNMAISIHPKKNSNEDIPPSTLLTYRWFLSYRMMTASMLCLCFSRQIN*"
  end
  it "should translate gene MhA1_Contig2992.frz3.gene1" do
    @contigsequence = @gff.sequencelist["MhA1_Contig2992"]
    @componentlist = {}
    @cdslist = {}
    @gff.each_CDS do | id, reclist, component |
      @componentlist[id] = component
      @cdslist[id] = reclist
    end
    name = "cds:MhA1_Contig2992.frz3.gene1"
    recs = @cdslist[name]
    component = @componentlist[name]
    # ntseq = @gff.assemble(@contigsequence,component.start,recs,:raw=>true)
    # ntseq.should == ""
    ntseq = @gff.assemble(@contigsequence,component.start,recs,:codonize=>true)
    ntseq.should == "AAAATTAATAAAAAAATAAATGATAATTCTTTTAATATTCAATCTGATTCGAATGAAAATTTGTTTAATGATGGAATTAATTCTGAACAAAATGAAGACAATATAGCAACAAAAAAAGGCAACAAAAAATTCGGTAAAAATCAAAAAGAAGGAAATAAAGAGTTGGATATTCAAAGTGAAGGTTTTGATAATAATGAAATACCTTCAAAAGAAAGCAAAAAACAAATAAGTAATTTTGGGGATAATGAAAGTGAATATGAAAAAGAAGAGGATAATAGAAAAAAGAAAGGGAAAAAAGGAATGATAGAAAAGTATGAATTAGGAAGGAATAAAGGAAGGGATAAAAATGAAAGAAATAAGGCTTCTGAAAGGTTTGATGAGCAGAATCAAGACAGAAATAATCAACGTGATAGTTTTGATTCTGGCAATAATGATAAATCACAAAGAGGCTTAGATAGCGGCACATTAGATGGAACAAATAATTTAAAAAGATCGAATGATGATCAATTACCAGAATTTTTGAAAACGGCCAGTCTCTCAGAGCGTCAGAAATTTCTTCAACTTGAAGCAGAAAATGACAGGTCCAAGTCTTCTATACGAAGAGATAAACAGAATTGGGCTGATCAACAAGGGCAGAGAATTTCTGATCTTTATAAACAATTTCAACAATCTTTACAACAAAAAGAAAAACAATTTAAAAGTGAACGTCAACGAAATGTTCAAATTAAATTAAGCAGAAATGCACAGAATGTTGATAAAAGAATTCAGGATCTTCTGAATAATCCTGATATTGCTGAAAGAGCTTTAATTCTTCAAATTGAACAAATCCTCGGCGGTACAGACGATAGTATTCGTCAGGAATTACAAAGACAAATATCTGTTATTGGACCATTAGATGGAAATATACCGCCAAATCTTACATAG"
    aaseq = @gff.assembleAA(@contigsequence,component.start,recs)
    aaseq.should == "KINKKINDNSFNIQSDSNENLFNDGINSEQNEDNIATKKGNKKFGKNQKEGNKELDIQSEGFDNNEIPSKESKKQISNFGDNESEYEKEEDNRKKKGKKGMIEKYELGRNKGRDKNERNKASERFDEQNQDRNNQRDSFDSGNNDKSQRGLDSGTLDGTNNLKRSNDDQLPEFLKTASLSERQKFLQLEAENDRSKSSIRRDKQNWADQQGQRISDLYKQFQQSLQQKEKQFKSERQRNVQIKLSRNAQNVDKRIQDLLNNPDIAERALILQIEQILGGTDDSIRQELQRQISVIGPLDGNIPPNLT*"
  end
  it "should fix Wormbase error MhA1_Contig3426.frz3.gene1" do
    @contigsequence = @gff.sequencelist["MhA1_Contig3426"]
    @componentlist = {}
    @cdslist = {}
    @gff.each_CDS do | id, reclist, component |
      @componentlist[id] = component
      @cdslist[id] = reclist
    end
    name = "cds:MhA1_Contig3426.frz3.gene1"
    recs = @cdslist[name]
    component = @componentlist[name]
    # :raw should not fix
    ntseq = @gff.assemble(@contigsequence,component.start,recs,:raw=>true)
    ntseq.should == "GCATCCAACAACAACAATTAGAAGTCTTTCCCAGCTCCTCCTCTGCCCCTCAGCAACAACAATACCCAGCGCAGCAGCTTCAATTAGTTACTCCTTTTATTGCATGCATAGCAGATGAATTGAGGGAGTTGATAGATGAAATGCGTATGTTTTAG"
    ntseq = @gff.assemble(@contigsequence,component.start,recs,:codonize=>true,:fix=>true)
    ntseq.should == "ATCCAACAACAACAATTAGAAGTCTTTCCCAGCTCCTCCTCTGCCCCTCAGCAACAACAATACCCAGCGCAGCAGCTTCAATTAGTTACTCCTTTTATTGCATGCATAGCAGATGAATTGAGGGAGTTGATAGATGAAATGCGTATGTTTTAG"
    ntseq.size.should == 153
    aaseq = @gff.assembleAA(@contigsequence,component.start,recs,:fix=>true)
    aaseq.should == "IQQQQLEVFPSSSSAPQQQQYPAQQLQLVTPFIACIADELRELIDEMRMF*"
  end
end



