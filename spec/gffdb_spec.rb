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

def iterators_should_be_implemented
  if TEST_NON_IMPLEMENTED
    it "should implement each_gene" 
    it "should implement each_gene" 
    it "should implement each_gene_seq" 
    it "should implement each_mRNA" 
    it "should implement each_exon" 
    it "should implement each_exon_seq" 
    it "should implement each_CDS" 
  end
  it "should implement each_mRNA_seq" do
    h = {} ; @gff.each_mRNA_seq { | id, seq | h[id] = seq }
    h["mrna01short Sequence:test01_1:400 (3:14)"].should == "GAAGATTTGTAT"
  end
  it "should implement each_CDS_seq" do
    h = {} ; @gff.each_CDS_seq { | id, seq | h[id] = seq }
    h.keys.sort[0].should == "cds1 Sequence:test01_1:400 (164:190, 192:200)"
    h["cds_short Sequence:test01_1:400 (3:14)"].should == "GAAGATTTGTAT"
  end
end

describe GFFdb, "GFF3 API with everything in memory" do

  before :all do 
    # initialize
    gffdb = Bio::GFFbrowser::GFFdb.new(TEST1)
    @gff = gffdb.assembler
  end

  iterators_should_be_implemented
end

describe GFFdb, "GFF3 API with :cache_components => 1000, :cache_records => :cache_none" do
  # iterators_should_be_implemented
end

describe GFFdb, "GFF3 API with :cache_components => 1000, :cache_records => 1000" do
  # iterators_should_be_implemented
end

describe GFFdb, "GFF3 API with :cache_records => :cache_none" do
  # iterators_should_be_implemented
end

describe GFFdb, "GFF3 API with :cache_components => :cache_none, :cache_records => :cache_none" do
  # iterators_should_be_implemented
end

describe GFFdb, "GFF3 API with external FASTA" do
  # iterators_should_be_implemented
end


