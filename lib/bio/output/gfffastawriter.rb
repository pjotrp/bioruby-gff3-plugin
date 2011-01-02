
module Bio
  module GFFbrowser

    class FastaWriter
      def initialize translate, validate
        @do_translate = translate
        @do_validate = validate
      end

      def put id, seq
        puts '>'+id
        put_seq seq
      end
      private

      def put_seq seq
        if @do_translate or @do_validate
          ntseq = Bio::Sequence::NA.new(seq)
          aaseq = ntseq.translate
          puts aaseq if @do_translate
          if @do_validate
            raise 'Validation problem of '+id if aaseq.count('*') > 1
          end
          return if @do_translate
        end
        puts seq
      end
    end
  end
end
