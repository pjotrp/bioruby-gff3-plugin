$: << '.'

require 'test/unit'
require 'regressiontest'

RegressionTest.create false

class Gff3Test < Test::Unit::TestCase

  rootpath = File.join(File.dirname(__FILE__),'..')
  p rootpath
  BIN = rootpath + '/bin/gff3-fetch'
  DAT = rootpath + '/test/data'
  def test_runs
    assert(RegressionTest.test `#{BIN} --logger stdout CDS --cache none #{DAT}/gff/test.gff3`,'test_gff3',"#{DAT}/regression")
    assert(RegressionTest.test `#{BIN} --logger stdout mRNA #{DAT}/gff/test-ext-fasta.fa #{DAT}/gff/test-ext-fasta.gff3`,'test_ext_gff3',"#{DAT}/regression")
  end
end
