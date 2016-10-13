package com.bina.varsim.types.variant;

//--- Java imports ---

import com.bina.intervalTree.SimpleInterval1D;
import com.bina.varsim.types.*;
import com.bina.varsim.util.SimpleReference;
import org.apache.log4j.Logger;

import java.io.UnsupportedEncodingException;
import java.util.Arrays;
import java.util.Random;

public class Variant implements Comparable<Variant>{
    private final static Logger log = Logger.getLogger(Variant.class.getName());
    public int idx = 0; // this is hopefully a unique index, for the split variants
    public int full_idx = 0; // this is hopefully a unique index, for the whole variants
    // this is the type before the variant was split into canonical ones
    public VariantOverallType original_type = null;
    // use a seed for reproducibility, should be an option or global
    private Random rand = null;
    private int pos = -1, referenceAlleleLength = -1;
    private byte[] ref;
    private ChrString chr;
    private FlexSeq[] alts;
    private byte maternal = 0, paternal = 0; // -1 for not avaliable
    private boolean isPhased = false; // Phasing
    private String filter;
    private String varId;
    // this is when the reference base is deleted
    // if it is the same as the first alt base
    //so refDeleted.length() <= 1 is always true?
    private String refDeleted;
    private String extraBase = "";
    private ChrString[] chr2;
    private int[] pos2;
    private int[] end2;
    private int end;
    private String[] translocationSubtype;

    public Variant(final Random rand) {
        // TODO define some methods to determine if a Variant is uninitialised
        this.rand = rand;
    }

    public Variant(final ChrString chr, final int pos, final int referenceAlleleLength, final byte[] ref,
                   final FlexSeq[] alts, final byte[] phase, final boolean isPhased, final String var_id, final String filter,
                   final String ref_deleted) {
        this(chr, pos, referenceAlleleLength, ref, alts, phase, isPhased, var_id, filter, ref_deleted, null, null, null, null, -1, null);
    }

    public Variant(ChrString chr, final int pos, final int referenceAlleleLength, final byte[] ref,
                   final FlexSeq[] alts, final byte[] phase, final boolean isPhased, final String var_id, final String filter,
                   final String ref_deleted, Random rand) {
        this(chr, pos, referenceAlleleLength, ref, alts, phase, isPhased, var_id, filter, ref_deleted, rand, null, null, null, -1, null);
    }
    public Variant(ChrString chr, final int pos, final int referenceAlleleLength, final byte[] ref,
                   final FlexSeq[] alts, final byte[] phase, final boolean isPhased, final String var_id, final String filter,
                   final String ref_deleted, Random rand, final ChrString[] chr2, final int[] pos2, final int[] end2, final int end, final String[] translocationSubtype) {
        this(rand);

        this.filter = filter;
        varId = var_id;

        this.chr = chr;
        this.pos = pos;
        this.referenceAlleleLength = referenceAlleleLength;

        // TODO we should put the reference matching code here
        this.ref = ref.clone();
        refDeleted = ref_deleted;
        this.alts = new FlexSeq[alts.length];
        for (int i = 0; i < alts.length; i++) {
            if (alts[i] != null) {
                this.alts[i] = new FlexSeq(alts[i]);
            } else {
                this.alts[i] = null;
            }
        }

        this.chr2 = chr2;
        this.pos2 = pos2;
        this.end2 = end2;
        paternal = phase[0];
        maternal = phase[1];
        this.isPhased = isPhased;
        this.end = end;
        this.translocationSubtype = translocationSubtype;
    }

    public Variant(final Variant var) {
        filter = var.filter;
        varId = var.varId;
        chr = var.chr;
        pos = var.pos;
        referenceAlleleLength = var.referenceAlleleLength;
        ref = var.ref.clone();
        refDeleted = var.refDeleted;
        alts = new FlexSeq[var.alts.length];
        for (int i = 0; i < var.alts.length; i++) {
            if (var.alts[i] != null) {
                alts[i] = new FlexSeq(var.alts[i]);
            } else {
                alts[i] = null;
            }
        }

        this.chr2 = var.chr2;
        this.pos2 = var.pos2;
        this.end2 = var.end2;
        paternal = var.paternal;
        maternal = var.maternal;
        isPhased = var.isPhased;
        rand = var.rand;
    }

    /**
     * @return Chromosome variant is on
     */
    public ChrString getChr() {
        if (chr == null) {
            throw new UnsupportedOperationException("ERROR: no legitimate chromosome name available!");
        }
        return chr;
    }

    /**
     * @return Start position of variant
     */
    public int getPos() {
        return pos;
    }

    // return false if fails
    // if return false, nothing is changed

    /**
     * Tries to move the variant to a new novel position. It checks if the new position is valid
     *
     * @param pos new novel position
     * @param ref reference sequence
     * @return return true if it was changed
     */
    public boolean setNovelPosition(final int pos, final SimpleReference ref) {

        // replace ref
        int len = this.ref.length;

        if (len > 0) {
            byte[] temp_ref = ref.byteRange(chr, pos, pos + len);

            for (byte b : temp_ref) {
                if (b == 'N') {
                    // don't allow N's
                    log.warn("N found at " + pos + " to " + (pos + len));
                    return false;
                }
            }

            for (FlexSeq f : alts) {
                if (f.getSeq() != null) {
                    // make sure there is no prefix the same
                    for (int i = 0; i < temp_ref.length; i++) {
                        if (i < f.getSeq().length) {
                            if (temp_ref[i] == f.getSeq()[i]) {
                                log.warn("Same ref at alt at " + pos + " to " + (pos + len));
                                return false;
                            } else {
                                break;
                            }
                        }
                    }

                }
            }

            this.ref = temp_ref;
        }

        // replace ref_deleted
        len = refDeleted.length();
        if (len > 0) {
            try {
                byte[] deleted_temp = ref.byteRange(chr, pos - len, pos);
                refDeleted = new String(deleted_temp, "US-ASCII");
            } catch (UnsupportedEncodingException e) {
                e.printStackTrace();
            }
        }

        // do the replacing
        this.pos = pos;

        return true;
    }

    /**
     * @param id variant id. usually the dbSNP id
     */
    public void setVarID(final String id) {
        varId = id;
    }

    /**
     * return the length of the getReferenceAlleleLength.
     * This is currently simply the length of the reference sequence that will
     * be replaced
     *
     * @return the length of the getReferenceAlleleLength
     */
    public int getReferenceAlleleLength() {
        return referenceAlleleLength;
    }

    /**
     * @return maternal allele index, 0 if reference
     */
    public int maternal() {
        return maternal;
    }

    /**
     * @return paternal allele index, 0 if reference
     */
    public int paternal() {
        return paternal;
    }

    /**
     * @param ind index of allele
     * @return the insertion sequence as a string
     */
    public byte[] insertion(final int ind) {
        if (ind <= 0 || ind > alts.length)
            return null;
        return alts[ind - 1].getSeq();
    }

    public ChrString getChr2(final int ind) {
        if (ind <= 0 || ind > alts.length)
            return new ChrString("");
        return this.chr2[ind - 1];
    }

    public int getPos2(final int ind) {
        if (ind <= 0 || ind > alts.length)
            return -1;
        return this.pos2[ind - 1];
    }

    public int getEnd2(final int ind) {
        if (ind <= 0 || ind > alts.length)
            return -1;
        return this.end2[ind - 1];
    }

    public ChrString[] getAllChr2() {
        return this.chr2;
    }

    public int[] getAllPos2() {
        return this.pos2;
    }

    public int[] getAllEnd2() {
        return this.end2;
    }

    public int getEnd() {
        return this.end;
    }
    public String[] getAllTranslocationSubtype() {
        return this.translocationSubtype;
    }

    /**
     * The length of an alternate allele, this is usually an insertion sequence
     * But in the case of SVs, not necessarily
     *
     * @param ind index of allele
     * @return the length of that allele
     */
    public int insertion_len(final int ind) {
        if (ind <= 0 || ind > alts.length)
            return 0;
        return alts[ind - 1].length();
    }

    // if it is a simple indel, it is just the length
    // if it is a complex variant, this is the maximum length of the insertion
    // and getReferenceAlleleLength
    public int maxLen(final int ind) {
        if (ind <= 0 || ind > alts.length)
            return 0;
        return Math.max(referenceAlleleLength, alts[ind - 1].length());
    }

    public int maxLen() {
        int len = 0;
        for (int i = 0; i < 2; i++) {
            len = Math.max(maxLen(get_allele(i)), len);
        }
        return len;
    }

    // this is the minimum length of the variants
    public int minLen() {
        int len = maxLen(get_allele(0));
        len = Math.min(maxLen(get_allele(1)), len);

        return len;
    }

    /*
    gets the interval enclosing the variant on the reference genome
    */
    public SimpleInterval1D get_interval(final int ind) {
        if (ind == 0 || referenceAlleleLength == 0) {
            return new SimpleInterval1D(pos, pos);
        }

        return new SimpleInterval1D(pos, pos + referenceAlleleLength - 1);
    }

    /*
    gets the interval for the variant, accounting for variant size
    */
    public SimpleInterval1D get_var_interval(final int ind) {
        try {
            if (ind == 0) {
                return new SimpleInterval1D(pos, pos);
            }

            // TODO hmm unsafe
            if (maxLen(ind) == Integer.MAX_VALUE) {
                return new SimpleInterval1D(pos, pos);
            } else {
                return new SimpleInterval1D(pos, pos + maxLen(ind) - 1);
            }
        } catch (RuntimeException e) {
            log.error("Bad variant interval: " + toString());
            log.error("pos: " + pos);
            log.error("ind: " + ind);
            log.error("maxLen(ind): " + maxLen(ind));
            e.printStackTrace();
            return null;
        }
    }

    // union of intervals from the genotypes
    public SimpleInterval1D get_geno_interval() {
        return get_interval(getgood_paternal()).union(get_interval(getgood_maternal()));
    }

    public SimpleInterval1D get_geno_var_interval() {
        return get_var_interval(getgood_paternal()).union(get_var_interval(getgood_maternal()));
    }

    public Genotypes getGeno() {
        return new Genotypes(getgood_paternal(), getgood_maternal());
    }

    /*
    * 0 = paternal
    * 1 = maternal
    * otherwise returns -1
     */
    public int get_allele(final int parent) {
        if (parent == 0) {
            return getgood_paternal();
        } else if (parent == 1) {
            return getgood_maternal();
        }
        return -1;
    }

    public void set_allele(final int parent, final byte allele) {
        if (parent == 0) {
            paternal = allele;
        } else if (parent == 1) {
            maternal = allele;
        }
    }

    // TODO this is wrong, but it only effects the count of variant bases
    public int variantBases() {
        int ret = referenceAlleleLength;
        for (FlexSeq _alt : alts) {
            if (referenceAlleleLength != _alt.length()) {
                ret += _alt.length();
            }
        }
        return ret;
    }

    /**
     * @param ind index of allele (starts at 1, 0 is reference)
     * @return type of allele at index ind
     */
    public VariantType getType(final int ind) {
        if (ind == 0) {
            return VariantType.Reference;
        }

        /*when variant type is explicitly specified in VCF,
        alt allele will be set to the correct type,
        we can be assured that correct variant type is
        returned.
         */
        FlexSeq.Type type = alts[ind - 1].getType();
        switch (type) {
            case DUP:
                return VariantType.Tandem_Duplication;
            case INS:
                return VariantType.Insertion;
            case INV:
                return VariantType.Inversion;
            case TRA:
                return VariantType.Translocation;
            default:
                break;
        }

        int inslen = insertion_len(ind);
        int dellen = referenceAlleleLength;
        if (inslen == 0 && dellen == 0) {
            return VariantType.Reference;
        } else if (inslen == 1 && dellen == 1) {
            return VariantType.SNP;
        } else if (inslen == 0 && dellen > 0) {
            return VariantType.Deletion;
        } else if (dellen - inslen > 0 && new String(ref).endsWith(alts[ind - 1].toString())) {
            return VariantType.Deletion;
        } else if (inslen > 0 && dellen == 0) {
            return VariantType.Insertion;
        } else if (inslen == dellen) {
            return VariantType.MNP;
        }
        return VariantType.Complex;
    }

    /**
     * @return overall type of the variant considering both alleles
     */
    public VariantOverallType getType() {
        int[] allele = {get_allele(0), get_allele(1)};

        // check Reference
        boolean is_ref = true;
        for (int a = 0; a < 2; a++) {
            if (getType(allele[a]) != VariantType.Reference) {
                is_ref = false;
                break;
            }
        }
        if (is_ref) {
            return VariantOverallType.Reference;
        }

        // check SNP
        boolean is_snp = true;
        for (int a = 0; a < 2; a++) {
            if (allele[a] > 0 && getType(allele[a]) != VariantType.SNP) {
                is_snp = false;
                break;
            }
        }
        if (is_snp) {
            return VariantOverallType.SNP;
        }

        // check INV
        boolean is_inv = true;
        for (int a = 0; a < 2; a++) {
            if (allele[a] > 0 && getType(allele[a]) != VariantType.Inversion) {
                is_inv = false;
                break;
            }
        }
        if (is_inv) {
            return VariantOverallType.Inversion;
        }

        // check DUP
        boolean is_dup = true;
        for (int a = 0; a < 2; a++) {
            if (allele[a] > 0 && getType(allele[a]) != VariantType.Tandem_Duplication) {
                is_dup = false;
                break;
            }
        }
        if (is_dup) {
            return VariantOverallType.Tandem_Duplication;
        }

        // check Deletion
        boolean is_del = true;
        for (int a = 0; a < 2; a++) {
            if (allele[a] > 0 && getType(allele[a]) != VariantType.Deletion) {
                is_del = false;
                break;
            }
        }
        if (is_del) {
            return VariantOverallType.Deletion;
        }

        // check DUP
        boolean is_ins = true;
        for (int a = 0; a < 2; a++) {
            if (allele[a] > 0 && getType(allele[a]) != VariantType.Insertion) {
                is_ins = false;
                break;
            }
        }
        if (is_ins) {
            return VariantOverallType.Insertion;
        }

        // check TRA
        boolean isTranslocation = true;
        for (int a = 0; a < 2; a++) {
            if (allele[a] > 0 && getType(allele[a]) != VariantType.Translocation) {
                isTranslocation = false;
                break;
            }
        }
        if (isTranslocation) {
            return VariantOverallType.Translocation;
        }


        /* Treat these as complex for now
        // check INDEL
        boolean is_indel = true;
        for (int a = 0; a < 2; a++) {
            if (allele[a] > 0 && !(getType(get_allele(a)) == Type.Deletion || getType(get_allele(a)) == Type.Insertion)) {
                is_indel = false;
                break;
            }
        }
        if (is_indel) {
            return VariantOverallType.INDEL;
        }

        // check MNP
        boolean is_mnp = true;
        for (int a = 0; a < 2; a++) {
            if (allele[a] > 0 && getType(get_allele(a)) != Type.MNP) {
                is_mnp = false;
                break;
            }
        }
        if (is_mnp) {
            return VariantOverallType.MNP;
        }
        */

        // otherwise it is complex
        return VariantOverallType.Complex;
    }

    /**
     * return alternative allele based on alternative allele index specified in GT field
     * alternative allele index = 1,2,...
     * @param ind
     * @return
     */
    public FlexSeq getAlt(final int ind) {
        if (ind <= 0 || ind > alts.length)
            return null;
        return alts[ind - 1];
    }

    public void setAlt(final int ind, FlexSeq alt) {
        if (ind <= 0 || ind > alts.length) {
            return;
        }
        alts[ind - 1] = alt;
    }

    public String getFilter() {
        return filter;
    }

    public boolean isPhased() {
        return isPhased;
    }

    public boolean isRef() {
        return paternal == 0 && maternal == 0;
    }

    public String getVar_id() {
        return varId;
    }

    public byte[] getRef() {
        return ref;
    }

    public String getOrig_Ref() {
        try {
            return refDeleted + new String(ref, "US-ASCII");
        } catch (UnsupportedEncodingException e) {
            e.printStackTrace();
            return "";
        }
    }

    public String getRef_deleted() {
        return refDeleted;
    }

    public int getCN(final int ind) {
        if (ind <= 0 || ind > alts.length)
            return 0;
        return alts[ind - 1].getCopy_num();
    }

    /**
     * @return true if any of the alternate alleles has copy number greater than 1
     */
    public boolean hasCN() {
        boolean CN_positive = false;
        for (FlexSeq _alt : alts) {
            if (_alt.getCopy_num() > 1) {
                CN_positive = true;
            }
        }

        return CN_positive;
    }

    public StringBuilder alt_string() {
        StringBuilder sbStr = new StringBuilder();
        for (int i = 0; i < alts.length; i++) {
            //if (i > 0 && alts[i].toString().equals(alts[i - 1].toString())) {
                /*Marghoob suggested that two identical symbolic alternative alleles are
                allowed. so essentially we go back to original behavior of VarSim.
                 */
            if (i > 0) {
                sbStr.append(",");
            }
            if (alts[i].isSeq()) {
                sbStr.append(refDeleted).append(alts[i].toString()).append(extraBase);
            } else {
                sbStr.append(alts[i].toString());
            }
        }
        return sbStr;
    }

    public int get_num_alt() {
        return alts.length;
    }

    /**
     * Randomly swap the haploype
     */
    public void randomizeHaplotype() {
        if (rand == null) {
            log.error("Cannot randomize haplotype");
            log.error(toString());
            System.exit(1);
        }

        if (rand.nextDouble() > 0.5) {
            return;
        }
        byte tmp = paternal;
        paternal = maternal;
        maternal = tmp;
    }

    /**
     * Randomize the genotype
     *
     * @param gender
     */
    public void randomizeGenotype(GenderType gender) {
        if (rand == null) {
            log.error("Cannot randomize genotype");
            log.error(toString());
            System.exit(1);
        }

        Genotypes g = new Genotypes(chr, gender, alts.length, rand);
        paternal = g.geno[0];
        maternal = g.geno[1];
    }


    /*
    Tests if all of the alternate alleles with sequence are ACTGN
     */
    public boolean isAltACTGN() {
        for (FlexSeq a : alts) {
            if (a.isSeq()) {
                if (!a.toString().matches("[ACTGN]*")) {
                    return false;
                }
            }
        }
        return true;
    }

    /*
    Returns true if the variant is homozygous
     */
    public boolean isHom() {
        return (paternal == maternal);
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        Variant variant = (Variant) o;

        if (referenceAlleleLength != variant.referenceAlleleLength) return false;
        if (isPhased != variant.isPhased) return false;
        if (maternal != variant.maternal) return false;
        if (paternal != variant.paternal) return false;
        if (pos != variant.pos) return false;
        if (full_idx != variant.full_idx) return false;
        if (idx != variant.idx) return false;
        if (!Arrays.equals(alts, variant.alts)) return false;
        if (chr != null ? !chr.equals(variant.chr) : variant.chr != null) return false;
        if (filter != null ? !filter.equals(variant.filter) : variant.filter != null) return false;
        if (rand != null ? !rand.equals(variant.rand) : variant.rand != null) return false;
        if (!Arrays.equals(ref, variant.ref)) return false;
        if (refDeleted != null ? !refDeleted.equals(variant.refDeleted) : variant.refDeleted != null)
            return false;
        if (varId != null ? !varId.equals(variant.varId) : variant.varId != null) return false;
        if (original_type != variant.original_type) return false;

        return true;
    }

    @Override
    public int hashCode() {
        int result = rand != null ? rand.hashCode() : 0;
        result = 31 * result + idx;
        result = 31 * result + full_idx;
        result = 31 * result + (original_type != null ? original_type.hashCode() : 0);
        result = 31 * result + pos;
        result = 31 * result + referenceAlleleLength;
        result = 31 * result + (ref != null ? Arrays.hashCode(ref) : 0);
        result = 31 * result + (chr != null ? chr.hashCode() : 0);
        result = 31 * result + (alts != null ? Arrays.hashCode(alts) : 0);
        result = 31 * result + (int) maternal;
        result = 31 * result + (int) paternal;
        result = 31 * result + (isPhased ? 1 : 0);
        result = 31 * result + (filter != null ? filter.hashCode() : 0);
        result = 31 * result + (varId != null ? varId.hashCode() : 0);
        result = 31 * result + (refDeleted != null ? refDeleted.hashCode() : 0);
        return result;
    }

    @Override
    public int compareTo(final Variant other) {
        final int chrCmp = chr.compareTo(other.chr);
        if (chrCmp != 0) {
            return chrCmp;
        }

        return getPos() - getRef_deleted().length() - (other.getPos() - other.getRef_deleted().length());
    }

    public String getLength() {
        StringBuilder len = new StringBuilder();

        for (int i = 0; i < alts.length; i++) {
            if (i > 0) {
                len.append(',');
            }

            VariantType t = getType(i + 1);

            if (t == VariantType.Deletion) {
                len.append(-referenceAlleleLength + alts[i].length()); // negative for deletions
            } else if (t == VariantType.Complex) {
                int alt_len = alts[i].length();
                if (referenceAlleleLength > alt_len) {
                    len.append(-referenceAlleleLength);
                } else {
                    len.append(alt_len);
                }
            } else {
                len.append(alts[i].length());
            }
        }

        return len.toString();
    }

    public void calculateExtraBase(final Sequence refSeq) {
        for (final FlexSeq alt : alts) {
            if (alt.isSeq() && alt.length() == 0 && getPos() + referenceAlleleLength < refSeq.length()) {
                //why extrabase is only 1-bp long?
                extraBase = String.valueOf((char) refSeq.byteAt(getPos() + referenceAlleleLength));
            }
        }
    }


    /**
     * @param sbStr will build a VCF record without genotype
     */
    // TODO, this should be self contained and output a VCF record
    private void buildVCFstr(final StringBuilder sbStr) {
        // chromosome name
        sbStr.append(chr.toString());
        sbStr.append("\t");
        // start position
        sbStr.append(pos - refDeleted.length());
        sbStr.append('\t');
        // variant id
        sbStr.append(varId);
        sbStr.append("\t");
        // ref allele
        sbStr.append(getOrig_Ref());
        sbStr.append("\t");
        // alt alleles
        sbStr.append(alt_string().toString());
        sbStr.append("\t");
        // variant quality
        sbStr.append(".\t");
        // pass label
        sbStr.append(filter);
        sbStr.append("\t");
        // INFO
        if (getType() == VariantOverallType.Tandem_Duplication) {
            sbStr.append("SVTYPE=DUP;");
            sbStr.append("SVLEN=");
            sbStr.append(getLength());
        } else if (getType() == VariantOverallType.Inversion) {
            sbStr.append("SVTYPE=INV;");
            sbStr.append("SVLEN=");
            sbStr.append(getLength());
        } else {
            sbStr.append("SVLEN=");
            sbStr.append(getLength());
        }
        sbStr.append("\t");

        // label (GT)
        if (hasCN()) {
            sbStr.append("CN:GT\t");
        } else {
            sbStr.append("GT\t");
        }

        if (hasCN()) {
            sbStr.append(String.valueOf(getCN(getgood_paternal())));
            sbStr.append("|");
            sbStr.append(String.valueOf(getCN(getgood_maternal())));
            sbStr.append(":");
        }

    }

    /**
     * @return a VCF record of the variant
     */
    public String toString() {
        StringBuilder sbStr = new StringBuilder();

        buildVCFstr(sbStr);


        sbStr.append(getgood_paternal());
        sbStr.append("|");
        sbStr.append(getgood_maternal());

        return sbStr.toString();
    }

    /**
     * @param paternal specified paternal allele
     * @param maternal specified maternal allele
     * @return the VCF record with prespecified genotype
     */
    public String toString(final int paternal, final int maternal) {
        StringBuilder sbStr = new StringBuilder();

        buildVCFstr(sbStr);

        // for this one we need to work out which one is added
        sbStr.append(paternal);
        sbStr.append("|");
        sbStr.append(maternal);

        return sbStr.toString();
    }

    public byte getgood_paternal() {
        if (paternal < 0) {
            return 1;
        } else {
            return paternal;
        }
    }

    public byte getgood_maternal() {
        if (maternal < 0) {
            return 1;
        } else {
            return maternal;
        }
    }

    public String getExtraBase() {
        return extraBase;
    }

}
