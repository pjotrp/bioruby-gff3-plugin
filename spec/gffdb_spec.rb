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

# GFF3 with included FASTA:
TESTGFF1='test/data/gff/test.gff3'
# GFF3 with external FASTA:
TESTGFF1EXT='test/data/gff/test-ext-fasta.gff3'
TESTGFF1FASTA='test/data/gff/test-ext-fasta.fa'

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

describe GFFdb, "GFF3 API (InMemory) with everything in memory" do

  before :all do 
    # initialize
    gffdb = Bio::GFFbrowser::GFFdb.new(TESTGFF1)
    @gff = gffdb.assembler
  end

  iterators_should_be_implemented
end

describe GFFdb, "GFF3 API with :cache_components => 1000, :cache_records => :cache_none" do
  # iterators_should_be_implemented
end

describe GFFdb, "GFF3 API with :cache_components => 1000, :cache_records => 1000" do
  it "should implement real caching"
  # iterators_should_be_implemented
end

describe GFFdb, "GFF3 API with :cache_records => :cache_none" do
  # iterators_should_be_implemented
end

describe GFFdb, "GFF3 API (NoCache) with :cache_components => :cache_none, :cache_records => :cache_none" do
  before :all do 
    # initialize
    gffdb = Bio::GFFbrowser::GFFdb.new(TESTGFF1, :cache_components => :cache_none, :cache_records => :cache_none)
    @gff = gffdb.assembler
  end

  iterators_should_be_implemented
end

describe GFFdb, "GFF3 API (InMemory) with external FASTA" do
  before :all do 
    # initialize
    gffdb = Bio::GFFbrowser::GFFdb.new(TESTGFF1EXT, :fasta_filename => TESTGFF1FASTA)
    @gff = gffdb.assembler
  end

  iterators_should_be_implemented
end

describe GFFdb, "GFF3 API (NoCache) with external FASTA" do
  it "should support external FASTA"
  # iterators_should_be_implemented
end


