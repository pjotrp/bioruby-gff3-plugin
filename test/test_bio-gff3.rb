$: << '.'
$:


require 'test/unit'

require "rubygems"
gem "bio-logger", ">= 0.6.0"
require 'bio-logger'

# require 'helper'
require 'regressiontest'

# Bio::Log::CLI.logger('stdout')
RegressionTest.create true


class Gff3Test < Test::Unit::TestCase

  def test_runs
    assert(RegressionTest.test `../bin/gff3-fetch --logger stdout CDS --cache none data/gff/test.gff3`,'test_gff3','data/regression')
    assert(RegressionTest.test `../bin/gff3-fetch --logger stdout mRNA data/gff/test-ext-fasta.fa data/gff/test-ext-fasta.gff3`,'test_ext_gff3','data/regression')
  end
end
