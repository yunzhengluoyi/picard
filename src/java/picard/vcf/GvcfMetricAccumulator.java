package picard.vcf;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import picard.util.DbSnpBitSetUtil;

import java.util.List;

/**
 * A accumulator for collecting metrics about a single-sample GVCF. The main point here is to sub-context
 * each {@link VariantContext} as it comes by to the alleles present in the genotype of the only sample.
 * Since this is a GVCF we expect a symbolic \<NON_REF\> allele to be present in each VC. This symbolic allele
 * will cause the regular {@link CallingMetricAccumulator} to return a small subset of the relevant metrics.
 */
public class GvcfMetricAccumulator extends CallingMetricAccumulator {
    String sample = null;

    public GvcfMetricAccumulator(final DbSnpBitSetUtil.DbSnpBitSets dbsnp) {
        super(dbsnp);
    }

    @Override
    public void setup(final VCFHeader vcfHeader) {
        final List<String> samples = vcfHeader.getGenotypeSamples();
        if (samples == null || samples.size() != 1) {
            throw new IllegalArgumentException("Expected to have exactly 1 sample in a GVCF, found " + ((samples == null) ? "0" : samples.size()));
        }
        sample = samples.get(0);
    }

    @Override
    public void accumulate(final VariantContext vc) {
        //since a gvcf always has a <NON_REF> allele, in order to get meaningful results we need to subcontext the variant to
        //the allele that actually appear in the only sample's genotype
        final VariantContext subContext = vc.subContextFromSample(sample);
        super.accumulate(subContext);
    }
}
