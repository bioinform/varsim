package com.binatechnologies.varsim;

//import org.apache.commons.lang3.StringUtils;
//import org.apache.log4j.Logger;

import java.util.ArrayList;

/**
 * Interval tree implemented as an AVL tree
 * This current implementation allows for duplicates
 * <p/>
 * Created by johnmu on 12/6/14.
 *
 * @author johnmu
 */
public class IntervalTree<Key extends Interval1D> {
    //private final static Logger log = Logger.getLogger(IntervalTree.class.getName());
    private IntervalTreeNode<Key> root = null;
    long numEntries = 0;

    public IntervalTree() {
    }

    /**
     * Initialises IntervalTree with some data
     *
     * @param data all of this will be added to the IntervalTree
     */
    public IntervalTree(Iterable<Key> data) {
        for (Key k : data) {
            add(k);
        }
    }

    /**
     * Initialise the tree with one node
     *
     * @param data data to add to the tree
     */
    public IntervalTree(Key data) {
        add(data);
    }

    /**
     * Some test code is here
     * TODO write actual unit tests
     *
     * @param args
     */
    public static void main(String args[]) {
        IntervalTree<Interval1D> test_tree = new IntervalTree<Interval1D>();
        System.out.println("Adding...");
        test_tree.add(new Interval1D(10, 20));
        test_tree.add(new Interval1D(5, 10));
        test_tree.add(new Interval1D(5, 20));
        test_tree.add(new Interval1D(15, 20));
        test_tree.testGetOverlap(new Interval1D(10, 10));
        test_tree.testGetOverlap(new Interval1D(11, 15));

        // test out the rotation
        // keep adding to the right and see if it balances
        test_tree.add(new Interval1D(15, 20));
        test_tree.add(new Interval1D(20, 20));
        test_tree.add(new Interval1D(20, 21));
        test_tree.add(new Interval1D(20, 24));
        test_tree.add(new Interval1D(30, 34));
        test_tree.add(new Interval1D(31, 34));
        test_tree.add(new Interval1D(35, 35));
        test_tree.add(new Interval1D(45, 45));
        test_tree.add(new Interval1D(55, 55));
        test_tree.add(new Interval1D(65, 65));
        test_tree.add(new Interval1D(75, 75));
        test_tree.add(new Interval1D(85, 85));

        // see what the tree looks like
        System.out.println("Print test tree...");
        test_tree.printTree();

        // keep adding to the left and see if it balances

        test_tree.add(new Interval1D(-10, -10));
        test_tree.add(new Interval1D(-20, -20));
        test_tree.add(new Interval1D(-30, -30));
        test_tree.add(new Interval1D(-40, -40));
        test_tree.add(new Interval1D(-50, -50));
        test_tree.add(new Interval1D(-60, -60));
        test_tree.add(new Interval1D(-70, -70));


        // see what the tree looks like
        System.out.println("Print test tree...");
        test_tree.printTree();

        test_tree.testGetOverlap(new Interval1D(20, 20));
        test_tree.testGetOverlap(new Interval1D(21, 22));
        test_tree.testGetOverlap(new Interval1D(30, 31));

        System.out.println("Do some rotation");
        System.out.println("**Rotate right: " + test_tree.rotateRight());
        System.out.println("**Rotate right: " + test_tree.rotateRight());
        System.out.println("**Rotate right: " + test_tree.rotateRight());
        System.out.println("**Rotate right: " + test_tree.rotateRight());
        System.out.println("**Rotate right: " + test_tree.rotateRight());
        System.out.println("**Rotate left: " + test_tree.rotateLeft());
        System.out.println("**Rotate left: " + test_tree.rotateLeft());
        System.out.println("**Rotate left: " + test_tree.rotateLeft());
        System.out.println("**Rotate left: " + test_tree.rotateLeft());
        System.out.println("**Rotate left: " + test_tree.rotateLeft());
        System.out.println("**Rotate left: " + test_tree.rotateLeft());
        System.out.println("**Rotate left: " + test_tree.rotateLeft());

        test_tree.testGetOverlap(new Interval1D(20, 20));
        test_tree.testGetOverlap(new Interval1D(21, 22));
        test_tree.testGetOverlap(new Interval1D(30, 31));

        test_tree.testContains(new Interval1D(20, 20));
        test_tree.testContains(new Interval1D(-45, -45));
    }

    private static void printList(final ArrayList<?> x) {
        for (Object val : x) {
            System.out.print(val + ",");
        }
        System.out.print('\n');
    }

    /**
     * Adds key to the root, ensures balance
     *
     * @param k key to add
     */
    public void add(Key k) {
        numEntries++;
        if (root == null) {
            root = new IntervalTreeNode<Key>(k);
        } else {
            add(root, null, k);
        }
    }

    /**
     * Adds key to a particular node in the tree
     *
     * @param head   Start adding at this node
     * @param parent null if head is the root, otherwise this is the parent of head
     * @param k      Key to add
     * @return The adjustment to balance factor
     */
    private int add(IntervalTreeNode<Key> head, IntervalTreeNode<Key> parent, Key k) {
        int balanceFactorChange = 0;
        if (head == null) {
            throw new RuntimeException("Tried to add to null node, key=" + k);
        } else {
            // Try to add the key
            int compVal = head.addKey(k);
            if (compVal == 0) {
                // added to center, we can stop now :D phew...
                return 0;
            } else if (compVal < 0) {
                if (head.getLeft() == null) {
                    head.setLeft(new IntervalTreeNode<Key>(k));
                    // Created a new branch, need to change the balance factor
                    balanceFactorChange--;
                } else {
                    balanceFactorChange = add(head.getLeft(), head, k);
                }
            } else {
                if (head.getRight() == null) {
                    head.setRight(new IntervalTreeNode<Key>(k));
                    // Created a new branch, need to change the balance factor
                    balanceFactorChange++;
                } else {
                    balanceFactorChange = add(head.getRight(), head, k);
                }
            }
        }
        // Change the balance factor as necessary
        if (balanceFactorChange > 0) {
            head.incBalanceFactor();
        } else if (balanceFactorChange < 0) {
            head.decBalanceFactor();
        }

        // Check the balance factor and rotate as necessary
        if (head.getBalanceFactor() > 1) {
            if (parent == null) {
                // we are at root
                root = rotateLeft(head);
            } else {
                parent.setChild(head, rotateLeft(head));
            }
            balanceFactorChange--;
        } else if (head.getBalanceFactor() < -1) {
            if (parent == null) {
                // we are at root
                root = rotateRight(head);
            } else {
                parent.setChild(head, rotateRight(head));
            }
            balanceFactorChange++;
        }

        return balanceFactorChange;
    }

    public ArrayList<Key> getOverlaps(Key k) {
        return getOverlaps(root, k, 0, 0);
    }

    public ArrayList<Key> getOverlaps(Key k, double reciprocalRatio) {
        return getOverlaps(root, k, reciprocalRatio, 0);
    }

    /**
     * Gets keys overlapping k from the root
     *
     * @param k Key to search for
     * @param reciprocalRatio
     * @param wiggle
     * @return ArrayList of keys that overlap the given Key
     */
    public ArrayList<Key> getOverlaps(Key k, double reciprocalRatio, int wiggle) {
        return getOverlaps(root, k, reciprocalRatio, wiggle);
    }

    /**
     * TODO There are some optimizations here if the left and right trees are augmented with more information
     *
     * @param k key to search for
     * @return List of Keys that overlap the given key
     */
    private ArrayList<Key> getOverlaps(IntervalTreeNode<Key> head, Key k, double reciprocalRatio, int wiggle) {
        if (head == null) {
            return new ArrayList<Key>();
        }
        ArrayList<Key> retVal = head.getOverlaps(k, reciprocalRatio, wiggle);
        int compVal = head.checkKey(k);
        if (compVal == 0) {
            retVal.addAll(getOverlaps(head.getLeft(), k, reciprocalRatio, wiggle));
            retVal.addAll(getOverlaps(head.getRight(), k, reciprocalRatio, wiggle));
        } else if (compVal < 0) {
            retVal.addAll(getOverlaps(head.getLeft(), k, reciprocalRatio, wiggle));
        } else {
            retVal.addAll(getOverlaps(head.getRight(), k, reciprocalRatio, wiggle));
        }
        return retVal;
    }

    public boolean contains(Key k) {
        return contains(root, k, 0, 0);
    }

    public boolean contains(Key k, double reciprocalRatio) {
        return contains(root, k, reciprocalRatio,0);
    }

    /**
     * @param k Key to search for
     * @param reciprocalRatio
     * @param wiggle
     * @return true if k overlaps and key in tree
     */
    public boolean contains(Key k, double reciprocalRatio, int wiggle) {
        return contains(root, k, reciprocalRatio, wiggle);
    }

    /**
     * @param head Start searching at this node
     * @param k    Key to search for
     * @return true if k overlaps and key from the tree rooted at head
     */
    private boolean contains(IntervalTreeNode<Key> head, Key k, double reciprocalRatio, int wiggle) {
        if (head == null) {
            return false;
        }
        if (head.contains(k, reciprocalRatio, wiggle)) {
            return true;
        }
        int compVal = head.checkKey(k);
        if (compVal == 0) {
            if (contains(head.getLeft(), k, reciprocalRatio, wiggle)) {
                return true;
            }
            if (contains(head.getRight(), k, reciprocalRatio, wiggle)) {
                return true;
            }
        } else if (compVal < 0) {
            if (contains(head.getLeft(), k, reciprocalRatio, wiggle)) {
                return true;
            }
        } else {
            if (contains(head.getRight(), k, reciprocalRatio, wiggle)) {
                return true;
            }
        }

        return false;
    }

    /**
     * rotate the root right and also adjust balance factor
     *
     * @return true if the root actually rotated
     */
    protected boolean rotateRight() {
        IntervalTreeNode<Key> prev_root = root;
        return prev_root != (root = rotateRight(root));
    }

    /**
     * Rotate the tree rooted at head right.
     * It will also adjust the balance factor
     *
     * @param head head of tree to rotate
     * @return new head of the tree, original head must be replaced or bad things will happen
     */
    private IntervalTreeNode<Key> rotateRight(IntervalTreeNode<Key> head) {
        if (head == null) {
            throw new RuntimeException("Tried to rotate null head right");
        }
        if (head.getLeft() == null) {
            // cannot rotate right in this case
            return head;
        }

        // compute the new balance factors
        int bfChild = head.getLeft().getBalanceFactor();
        int bfHead = head.getBalanceFactor();
        int bfNewHead;
        int bfNewChild;

        if (bfChild >= 0) {
            bfNewChild = bfHead + 1;
            if (bfNewChild >= 0) {
                bfNewHead = bfHead + bfChild + 2;
            } else {
                bfNewHead = bfChild + 1;
            }
        } else {
            bfNewChild = bfHead - bfChild + 1;
            if (bfNewChild >= 0) {
                bfNewHead = bfHead + 2;
            } else {
                bfNewHead = bfChild + 1;
            }
        }

        IntervalTreeNode<Key> newHead = head.getLeft();
        newHead.setBalanceFactor(bfNewHead);
        head.setBalanceFactor(bfNewChild);
        IntervalTreeNode<Key> temp = newHead.getRight();
        newHead.setRight(head);
        head.setLeft(temp);
        return newHead;
    }

    /**
     * Rotate the root left and also adjust the balance factor
     *
     * @return true if the root actually rotated
     */
    protected boolean rotateLeft() {
        IntervalTreeNode<Key> prev_root = root;
        return prev_root != (root = rotateLeft(root));
    }

    /**
     * Rotate the tree rooted at head left.
     * It will also adjust the balance factor
     *
     * @param head head of tree to rotate
     * @return new head of the tree, original head must be replaced or bad things will happen
     */
    private IntervalTreeNode<Key> rotateLeft(IntervalTreeNode<Key> head) {
        if (head == null) {
            throw new RuntimeException("Tried to rotate null head left");
        }
        if (head.getRight() == null) {
            // cannot rotate left
            return head;
        }

        // compute the new balance factors
        int bfChild = head.getRight().getBalanceFactor();
        int bfHead = head.getBalanceFactor();
        int bfNewHead;
        int bfNewChild;

        if (bfChild <= 0) {
            // left branch on child is heavier
            bfNewChild = bfHead - 1;
            if (bfNewChild >= 0) {
                // left branch of child is heavier than left branch of head
                bfNewHead = bfChild - 1;
            } else {
                bfNewHead = bfHead + bfChild - 2;
            }
        } else {
            bfNewChild = bfHead - bfChild - 1;
            if (bfNewChild >= 0) {
                // left branch of child is heavier than left branch of head
                bfNewHead = bfChild - 1;
            } else {
                bfNewHead = bfHead - 2;
            }
        }

        IntervalTreeNode<Key> newHead = head.getRight();
        newHead.setBalanceFactor(bfNewHead);
        head.setBalanceFactor(bfNewChild);
        IntervalTreeNode<Key> temp = newHead.getLeft();
        newHead.setLeft(head);
        head.setRight(temp);
        return newHead;
    }

    public long size(){
        return numEntries;
    }

    /**
     * Ugly print the whole tree
     */
    public void printTree() {
        printTree(root, 0);
    }

    /**
     * Ugly print the tree rooted at head
     *
     * @param head  start printing from this node
     * @param level level head is at, this is just for printing
     */
    private void printTree(final IntervalTreeNode<Key> head, int level) {
        System.out.println(level + ":" + head);
        if (head != null) {
            System.out.println("Right: " + (level + 1));
            printTree(head.getRight(), level + 1);
            System.out.println("Left: " + (level + 1));
            printTree(head.getLeft(), level + 1);
        }
    }

    private void testGetOverlap(Key k) {
        System.out.println("Getting..." + k);
        //System.out.println(StringUtils.join(output, ','));
        printList(getOverlaps(k));
    }

    private void testContains(Key k) {
        System.out.println("Contains..." + k + " -- " + contains(k));
    }
}
