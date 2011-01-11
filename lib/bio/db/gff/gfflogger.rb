
module Bio
  module GFFbrowser

    module Helpers

      module Error
        include Bio::Log

        def info str, id=''
          log = LoggerPlus['bio-gff3']
          log.info str+" <#{id}>"
        end

        def warn str, id=''
          log = LoggerPlus['bio-gff3']
          log.warn str+" <#{id}>"
        end

        def error str, id=''
          log = LoggerPlus['bio-gff3']
          log.error str+" <#{id}>"
        end
      end
    end
  end
end
