# RSpec for BioRuby-GFF3-Plugin. Run with something like:
#
#   ruby -I ../bioruby/lib/ ~/.gems/bin/spec spec/gff3_spec.rb 
#
# Copyright (C) 2010 Pjotr Prins <pjotr.prins@thebird.nl>
#
$: << "../lib"

require 'bio/db/gff/gffdb'

TEST1='test/data/gff/test.gff3'
TEST2='test/data/gff/standard.gff3'

describe Bio::GFF::GFF3::FileIterator, "iterates a file" do

  before :all do 
    # initialize
    @iter = Bio::GFF::GFF3::FileIterator.new(TEST1)
  end

  it "should parse a file and yield records" do 
    @iter.each_rec do | id, rec |
      # p [id, rec, rec.io_seek]
      rec.io_seek.should == 51
      break
    end
  end

  it "should handle embedded FASTA records" do
    @iter = Bio::GFF::GFF3::FileIterator.new(TEST1)
    last = nil
    @iter.each_rec(Proc.new { | rec | Bio::GFFbrowser::Helpers::Record::formatID(rec) }) do | id, rec |
      p [id, rec]
      last = rec
    end
    last.source[0..5].should == "test01"
    firstid = 'unknown'
    @iter.each_sequence do | id, seq |
      firstid = id
    end
    firstid.should == "test01"
  end
end


