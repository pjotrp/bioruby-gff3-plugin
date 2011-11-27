# Regression tester
#
# Info:: Pjotr's shared Ruby modules
# Author:: Pjotr Prins
# mail:: pjotr.prins@thebird.nl
# Copyright:: July 2007
# License:: Ruby License

module RegressionTest

	def RegressionTest.create b
		@@test_create = b
	end
	
	#  Invoke the regression test by passing a string - which ends up a file
	#  in test/regression with +filename+. When +create+ is +true+ the file
	#  will be created/overwritten. Otherwise it is tested against returning
	#  whether it has equal or not. When a test fails both test file and new
	#  file exist in the regression directory - so you can execute a diff.
	#
	#  Example:
	#    RegressionTest.test `#{cfrubybin} --help`,'cfruby_helptext',$test_create

	def RegressionTest.test text, filename, testdir, create = @@test_create
		fn = testdir+'/'+filename+'.rtest'
		fntest = fn+'.new'
		
		if create
			f = File.open(fn,'w')
			f.write text
			File.unlink fntest if File.exist? fntest
		else
			# ---- here we have to compare info
			if ! File.exist?(fn)
				raise "Cannot execute regression test because file #{fn} does not exist! - use --create option?"
			end
			f = File.open(fn)
			b = ''
			f.each do | line |
				b += line
			end
			if b!=text
				# ---- Write newer file
				f2 = File.open(fntest,'w')
				f2.write text
				return false
			end
		end
		true
	end

end
