# find local plugin installation, and use it when there
rootpath = File.dirname(File.dirname(__FILE__))
bio_logger_path = File.join(rootpath,'..','bioruby-logger','lib')
if File.directory? bio_logger_path
  $: << bio_logger_path
	$stderr.print "bio-logger loaded directly\n"
else
  require "rubygems"
  gem "bio-logger", ">= 0.6.0"
end
require 'bio-logger'

Bio::Log::LoggerPlus.new('bio-gff3')

require 'bio'
require 'bio/db/gff/gfflogger'
require 'bio/db/gff/gffvalidate'
require 'bio/db/gff/gffsection'
require 'bio/db/gff/gffcomponent'
require 'bio/db/gff/gffsequence'
require 'bio/db/gff/gfffileiterator'
require 'bio/db/gff/gfffasta'
require 'bio/db/gff/gffparser'
require 'bio/db/gff/gffinmemory'
require 'bio/db/gff/gffnocache'
require 'bio/db/gff/gffdb'

