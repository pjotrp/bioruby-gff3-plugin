# RSpec for BioRuby-GFF3-Plugin. Run with something like:
#
#   ruby -I ../bioruby/lib/ ~/.gems/bin/spec spec/gffdb_spec.rb 
#
# Copyright (C) 2010 Pjotr Prins <pjotr.prins@thebird.nl>
#
$: << "../lib"


require 'bio/db/gff/gffdb'

include Bio::GFFbrowser

TEST_NON_IMPLEMENTED=false

TEST1='test/data/gff/test.gff3'

describe GFFdb, "GFF3 API with :cache_none" do

  before :all do 
    # initialize
    gff3 = Bio::GFF::GFF3.new(File.read(TEST1))
    @gffdb = Bio::GFFbrowser::GFFdb.new(gff3)
  end

  if TEST_NON_IMPLEMENTED
    it "should implement each_gene" 
    it "should implement each_gene" 
    it "should implement each_gene_seq" 
    it "should implement each_mRNA" 
    it "should implement each_exon" 
    it "should implement each_exon_seq" 
    it "should implement each_CDS" 
  end
  it "should implement each_mRNA_seq" 
  # do
  #   seq.should == "GAAGATTTGTAT"
  # end


  it "should implement each_CDS_seq" do
    res = []
    @gffdb.each_CDS_seq do | id, seq |
      res << id
    end
    res[0].should == "cds2 Sequence:test01_1:400 (192:200)"
  end
end

describe GFFdb, "GFF3 API with :cache_all" do
  if TEST_NON_IMPLEMENTED
    it "should implement each_gene" 
    it "should implement each_gene" 
    it "should implement each_gene_seq" 
    it "should implement each_mRNA" 
    it "should implement each_mRNA_seq" 
    it "should implement each_exon" 
    it "should implement each_exon_seq" 
    it "should implement each_CDS" 
    it "should implement each_CDS_seq" 
  end
end

describe GFFdb, "GFF3 API with :cache_component" do
  if TEST_NON_IMPLEMENTED
    it "should implement each_gene" 
    it "should implement each_gene" 
    it "should implement each_gene_seq" 
    it "should implement each_mRNA" 
    it "should implement each_mRNA_seq" 
    it "should implement each_exon" 
    it "should implement each_exon_seq" 
    it "should implement each_CDS" 
    it "should implement each_CDS_seq" 
  end
end

describe GFFdb, "GFF3 API with external FASTA" do
  if TEST_NON_IMPLEMENTED
    it "should implement each_gene" 
    it "should implement each_gene" 
    it "should implement each_gene_seq" 
    it "should implement each_mRNA" 
    it "should implement each_mRNA_seq" 
    it "should implement each_exon" 
    it "should implement each_exon_seq" 
    it "should implement each_CDS" 
    it "should implement each_CDS_seq" 
  end
end


