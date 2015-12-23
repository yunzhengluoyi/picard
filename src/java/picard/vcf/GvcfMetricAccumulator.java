package picard.vcf;

import htsjdk.variant.variantcontext.VariantContext;
import picard.util.DbSnpBitSetUtil;

import java.util.List;

/**
 * Created by farjoun on 12/21/15.
 */
public class GvcfMetricAccumulator extends CallingMetricAccumulator {
    String sample = null;
    public GvcfMetricAccumulator(final DbSnpBitSetUtil.DbSnpBitSets dbsnp) {
        super(dbsnp);
    }


    @Override
    public void accumulate(final VariantContext vc) {
        if (sample == null){
            final List<String> samples = vc.getSampleNamesOrderedByName();
            if (samples == null || samples.size() != 1){
                throw new IllegalArgumentException("Expected to have exactly 1 sample in a GVCF, found " + ((samples == null) ? "NULL" : samples.size()));
            }
            sample = samples.get(0);
        }
        VariantContext subContext = vc.subContextFromSample(sample);
        super.accumulate(subContext);
    }
}
