# RSpec for BioRuby-GFF3-Plugin. Run with something like:
#
#   ruby -I ../bioruby/lib/ ~/.gems/bin/spec spec/gff3_spec.rb 
#
# Copyright (C) 2010 Pjotr Prins <pjotr.prins@thebird.nl>
#
$: << "../lib"

require 'bio/db/gff/gffdb'

TEST1='test/data/gff/test.gff3'

describe Bio::GFF::GFF3::FileIterator, "iterates a file" do

  before :all do 
    # initialize
    iter = Bio::GFF::GFF3::FileIterator.new(TEST1)
  end

  it "should parse a file and yield records"
end


