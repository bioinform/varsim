package com.bina.varsim.fastqLiftover.types;

public class MapBlock implements Comparable<MapBlock> {
    public int size;
    public GenomeLocation srcLoc;
    public GenomeLocation dstLoc;
    public BlockType blockType;
    public int direction;
    public String name;

    public MapBlock(final int size, final String srcChr, final int srcLocation, final String dstChr, final int dstLocation, final String direction, final String featureType, final String name) {
        this.size = size;
        srcLoc = new GenomeLocation(srcChr, srcLocation);
        dstLoc = new GenomeLocation(dstChr, dstLocation);
        this.direction = direction.equals("+") ? 0 : 1;
        this.blockType = BlockType.fromName(featureType);
        this.name = name;
    }

    public MapBlock(final GenomeLocation srcLoc) {
        size = 0;
        this.srcLoc = srcLoc;
        this.dstLoc = new GenomeLocation("", 0);
        this.blockType = BlockType.UNKNOWN;
        this.direction = 0;
    }

    @Override
    public int compareTo(MapBlock rhs) {
        return srcLoc.location - rhs.srcLoc.location;
    }

    public String toString() {
        return srcLoc.toString() + " -> " + dstLoc.toString() + " " + blockType.name() + ((direction == 0) ? "+" : "-") + " " + size + "bp";
    }

    @Override
    public boolean equals(Object rhs) {
        boolean result = false;
        if (rhs instanceof MapBlock) {
            MapBlock rhsMapBlock = (MapBlock) rhs;
            result = (size == rhsMapBlock.size && srcLoc.chromosome.equals(rhsMapBlock.srcLoc.chromosome) && srcLoc.location == rhsMapBlock.srcLoc.location);
        }
        return result;
    }

    public enum BlockType {
        SEQ("S", "Sequence"), INS("I", "Insertion"), DEL("D", "Deletion"), INV("V", "Inversion"), DUP_TANDEM("T", "Tandem_Duplication"), UNKNOWN("U", "Unknown");

        private final String shortName;
        private final String longName;

        BlockType(String shortName, String longName) {
            this.shortName = shortName;
            this.longName = longName;
        }

        public static BlockType fromName(final String s) {
            for (final BlockType blockType : values()) {
                if (blockType.name().equals(s) || blockType.shortName.equals(s) || blockType.longName.equals(s)) {
                    return blockType;
                }
            }
            return "".equals(s) ? SEQ : UNKNOWN;
        }

        public boolean equalsName(String otherName) {
            return (otherName == null) ? false : shortName.equals(otherName);
        }

        public String toString() {
            return shortName;
        }

        public String getShortName() {
            return shortName;
        }

        public String getLongName() {
            return longName;
        }
    }
}
