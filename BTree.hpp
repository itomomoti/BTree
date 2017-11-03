/*!
 * Copyright (c) 2017 Tomohiro I
 *
 * This program is released under the MIT License.
 * http://opensource.org/licenses/mit-license.php
 */
/*!
 * @file BTree.hpp
 * @brief Pointer-based implementation of upper part of B+tree.
 * @attention Bottom part of B+tree can be implemented in more space efficient way (e.g., array-based implementation).
 * @author Tomohiro I
 * @date 2017-01-26
 * @todo Support delete.
 */
#ifndef INCLUDE_GUARD_BTree
#define INCLUDE_GUARD_BTree

#include <algorithm>
#include <cassert>

namespace itmmti
{
  /*!
   * @brief Pointer-based implementation of upper part of B+tree.
   * @attention Bottom part of B+tree can be implemented in more space efficient way (e.g., array-based implementation).
   * @tparam kB should be in {4, 8, 16, 32, 64, 128}. kB/2 <= 'numChildren_' <= kB.
   * @tparam ROW_NUM: We can maintain 'ROW_NUM' different types of psum of weights.
   */
  template <uint8_t kB = 64, uint8_t ROW_NUM = 1>
  class BTreeNode
  {
  private:
    // Private constant, alias etc.
    using BTreeNodeT = BTreeNode<kB, ROW_NUM>;

    /*!
     * @brief For representing current status of node by 'flags_'.
     */
    enum {
      isBorderBit = 1, //!< This node lies in border, i.e., its children are bottom nodes.
      isRootBit = 2, //!< This node is a root.
      isDummyBit = 4, //!< This node is a dummy node.
    };


  public:
    // Public constant, alias etc.
    enum {
      NOTFOUND = UINTPTR_MAX
    };


  private:
    // Member variables
    uint64_t psum_[ROW_NUM][kB+1]; //!< Partial sum: psum_[i+1] = sum_{i = 0}^{i} [weight of i-th child (0base)]
    BTreeNodeT * parent_; //!< Pointer to parent node.
    uint8_t idxInSibling_; //!< This node is 'idxInSibling_'-th child (0base) of its parent.
    uint8_t numChildren_; //!< Current num of children.
    uint8_t flags_;
    BTreeNodeT * children_[kB]; //!< Pointers to children. Note that the children might not BTreeNode type when this node is in border.
    BTreeNodeT * lmBorderNode_; //!< To remember leftmost border node of this node.


  public:
    BTreeNode<kB, ROW_NUM>(bool isBorder, bool isRoot, BTreeNode<kB, ROW_NUM> * lmBorderNode, bool isDummy = false)
    : parent_(NULL),
      numChildren_(0),
      flags_(isBorder * isBorderBit | isRoot * isRootBit | isDummy * isDummyBit),
      lmBorderNode_(lmBorderNode)
    {
      psum_[0] = 0;
    }
    ~BTreeNode<kB, ROW_NUM>() = default;
    BTreeNode<kB, ROW_NUM>(const BTreeNode<kB, ROW_NUM> &) = delete;
    BTreeNode<kB, ROW_NUM> & operator=(const BTreeNode<kB, ROW_NUM> &) = delete;


    /*!
     * @brief Clear all upper nodes in B+tree.
     * @node Bottom nodes might be cleared manually before calling this function.
     */
    void clearUpperNodes() noexcept {
      if (!isBorder()) {
        for (uint8_t i = 0; i < numChildren_; ++i) {
          children_[i]->clearTree();
        }
      }
      delete this;
    }


    //// simple getter
  public:
    /*!
     * @brief Get the sum of weights (of row number "ROW") of child subtrees numbered from 0 to i-1 (note that i is NOT included).
     */
    template<uint8_t ROW>
    uint64_t getPSum
    (
     const uint8_t idx_excl //!< in [0.."numChildren_"].
     ) const noexcept {
      static_assert(ROW < ROW_NUM, "ROW should be smaller than ROW_NUM");
      assert(idx_excl <= numChildren_);

      return psum_[ROW][idx_excl];
    }


    /*!
     * @brief Get the weight (of row number "ROW") of a child subtree.
     */
    template<uint8_t ROW>
    uint64_t getWeightOfChild
    (
     const uint8_t idx //!< in [0.."numChildren_").
     ) const noexcept {
      static_assert(ROW < ROW_NUM, "ROW should be smaller than ROW_NUM");
      assert(idx < numChildren_);

      return psum_[ROW][idx+1] - psum_[ROW][idx];
    }


    /*!
     * @brief Get the weights of a child subtree.
     */
    void getWeightOfChild
    (
     const uint8_t idx, //!< in [0.."numChildren_").
     uint64_t (&weights)[ROW_NUM] //!< [out]
     ) const noexcept {
      assert(idx < numChildren_);

      for (uint8_t row = 0; row < ROW_NUM; ++row) {
        weights[row] = psum_[row][idx+1] - psum_[row][idx];
      }
    }


    /*!
     * @brief Get const psum array.
     */
    template<uint8_t ROW>
    const uint64_t * getConstPSumArray () const noexcept {
      return psum_[ROW];
    }


    /*!
     * @brief Get the weight (of row number "row") of this node.
     */
    template<uint8_t ROW>
    uint64_t getSumOfWeight() const noexcept {
      static_assert(ROW < ROW_NUM, "ROW should be smaller than ROW_NUM");

      return psum_[ROW][numChildren_];
    }


    /*!
     * @brief Get pointer to the i-th child (0base).
     */
    BTreeNodeT * getChildPtr
    (
     uint8_t i //!< in [0, 'numChildren_').
     ) const noexcept {
      assert(i < numChildren_);

      return children_[i];
    }


    /*!
     * @brief Get pointer to parent.
     */
    BTreeNodeT * getParent() const noexcept {
      return parent_;
    }


    /*!
     * @brief Get BTreeNode<kB, ROW_NUM>::idxInSibling_.
     */
    uint8_t getIdxInSibling() const noexcept {
      return idxInSibling_;
    }


    /*!
     * @brief Get num of children.
     */
    uint8_t getNumChildren() const noexcept {
      return numChildren_;
    }


    /*!
     * @brief Return if this node is in border.
     */
    bool isBorder() const noexcept {
      return flags_ & isBorderBit;
    }


    /*!
     * @brief Return if this node is root.
     */
    bool isRoot() const noexcept {
      return flags_ & isRootBit;
    }


    /*!
     * @brief Return if this node is a dummy node.
     */
    bool isDummy() const noexcept {
      return flags_ & isDummyBit;
    }


  public:
    //// function to traverse tree
    /*!
     * @brief Get leftmost border node. It takes O(1) time.
     */
    BTreeNodeT * getLmBorderNode() const noexcept {
      return lmBorderNode_;
    }


    /*!
     * @brief Get rightmost border node. It takes O(h) time, where h is the height of this node.
     */
    BTreeNodeT * getRmBorderNode() const noexcept {
      const auto * node = this;
      while (!(node->isBorder())) {
        node = node->children_[node->getNumChildren() - 1];
      }
      return node;
    }


    /*!
     * @brief Get leftmost bottom. It takes O(1) time.
     */
    BTreeNodeT * getLmBtm() const noexcept {
      return lmBorderNode_->getChildPtr(0);
    }


    /*!
     * @brief Get rightmost bottom. It takes O(h) time, where h is the height of this node.
     */
    BTreeNodeT * getRmBtm() const noexcept {
      const auto * node = this;
      while (!(node->isBorder())) {
        node = node->children_[node->getNumChildren() - 1];
      }
      return node->children_[node->getNumChildren() - 1];
    }


    /*!
     * @brief Return next bottm starting from "idxInSib"-th child (0base) of this node.
     */
    BTreeNodeT * getNextBtm
    (
     uint8_t idxInSib
     ) const noexcept {
      const auto * node = this;
      while (idxInSib + 1 == node->getNumChildren() && !(node->isRoot())) {
        idxInSib = node->getIdxInSibling();
        node = node->getParent();
      }
      if (idxInSib + 1 < node->getNumChildren()) {
        if (node->isBorder()) {
          return node->getChildPtr(idxInSib + 1);
        } else {
          return node->getChildPtr(idxInSib + 1)->getLmBtm();
        }
      }
      return reinterpret_cast<BTreeNodeT *>(NOTFOUND);
    }


    /*!
     * @brief Return previous bottm starting from "idxInSib"-th child (0base) of this node.
     */
    BTreeNodeT * getPrevBtm
    (
     uint8_t idxInSib
     ) const noexcept {
      const auto * node = this;
      while (idxInSib == 0 && !(node->isRoot())) {
        idxInSib = node->getIdxInSibling();
        node = node->getParent();
      }
      if (idxInSib) {
        if (node->isBorder()) {
          return node->getChildPtr(idxInSib - 1);
        } else {
          return node->getChildPtr(idxInSib - 1)->getRmBtm();
        }
      }
      return reinterpret_cast<BTreeNodeT *>(NOTFOUND);
    }


    /*!
     * @brief Return next bottm starting from "idxInSib"-th child (0base) of this node.
     */
    BTreeNodeT * getNextBtmRef
    (
     uint8_t & idxInSib //!< [in,out]
     ) const noexcept {
      const auto * node = this;
      while (++idxInSib == node->getNumChildren() && !(node->isRoot())) {
        idxInSib = node->getIdxInSibling();
        node = node->getParent();
      }
      if (idxInSib < node->getNumChildren()) {
        return node;
      }
      return reinterpret_cast<BTreeNodeT *>(NOTFOUND);
    }


    /*!
     * @brief Return previous bottm starting from "idxInSib"-th child (0base) of this node.
     */
    BTreeNodeT * getPrevBtmRef
    (
     uint8_t & idxInSib //!< [in,out]
     ) const noexcept {
      const auto * node = this;
      while (idxInSib == 0 && !(node->isRoot())) {
        idxInSib = node->getIdxInSibling();
        node = node->getParent();
      }
      if (idxInSib) {
        --idxInSib;
        return node;
      }
      return reinterpret_cast<BTreeNodeT *>(NOTFOUND);
    }


    /*!
     * @brief
     *   Return partial sum (or row number "ROW") up to the node
     *   (inclusive iff "inclusive == true") indicated by "idx"-th child (0base) of this node.
     */
    /*
      template<uint8_t ROW>
      uint64_t calcPSum
      (
      uint8_t idx,
      bool inclusive = false
      ) const noexcept {
      static_assert(ROW < ROW_NUM, "ROW should be smaller than ROW_NUM");
      assert(isBorder());

      const auto * node = this;
      uint64_t ret = 0;
      while (true) {
      ret += node->getPSum<ROW>(idx + inclusive);
      if (node->isRoot()) {
      return ret;
      }
      idx = node->getIdxInSibling();
      node = node->getParent();
      }
      }
    */


    /*!
     * @brief Return partial sum (of row number "ROW") up to the node (exclusive) indicated by "idx"-th child (0base) of this node.
     */
    template<uint8_t ROW>
    uint64_t calcPSum
    (
     uint8_t idx
     ) const noexcept {
      static_assert(ROW < ROW_NUM, "ROW should be smaller than ROW_NUM");
      assert(isBorder());

      const auto * node = this;
      uint64_t ret = 0;
      while (true) {
        ret += node->getPSum<ROW>(idx);
        if (node->isRoot()) {
          return ret;
        }
        idx = node->getIdxInSibling();
        node = node->getParent();
      }
    }


    /*!
     * @brief Return partial sum (of row number "ROW") up to the node (exclusive) indicated by "idx"-th child (0base) of this node.
     */
    template<uint8_t ROW>
    uint64_t calcPSum
    (
     uint8_t idx,
     BTreeNodeT * retNode //!< [out] To capture the root of the BTree.
     ) const noexcept {
      static_assert(ROW < ROW_NUM, "ROW should be smaller than ROW_NUM");
      assert(isBorder());

      retNode = this;
      uint64_t ret = 0;
      while (true) {
        ret += retNode->getPSum<ROW>(idx);
        if (retNode->isRoot()) {
          return ret;
        }
        idx = retNode->getIdxInSibling();
        retNode = retNode->getParent();
      }
    }


    /*!
     * @brief Traverse tree looking for "pos" in weights array of "ROW".
     * @return Pointer to bottom node where partial sum of "pos" is achieved, where weight-0 nodes (e.g. dummy nodes) are skipped.
     */
    template<uint8_t ROW>
    BTreeNodeT * searchPos
    (
     uint64_t & pos //!< [in,out] Give global position to search. It is modified to relative position in bottom node.
     ) const noexcept {
      static_assert(ROW < ROW_NUM, "ROW should be smaller than ROW_NUM");
      assert(pos < this->getSumOfWeight<ROW>());

      BTreeNodeT * node = this;
      while (true) {
        uint8_t i = 0;
        auto array = node->getConstPSumArray<ROW>();
        while (pos >= array[i + 1]) {
          ++i;
        }
        pos -= array[i];
        if (isBorder()) {
          return children_[i];
        }
        node = node->getChildPtr(i);
      }
    }


    /*!
     * @brief Traverse tree looking for "pos" in weights array of "ROW".
     * @return Pointer to bottom node where partial sum of "pos" is achieved, where weight-0 nodes (e.g. dummy nodes) are skipped.
     */
    template <uint8_t ROW>
    void searchPos
    (
     uint64_t & pos, //!< [in,out] Give global position to search. It is modified to relative position in bottom node.
     BTreeNodeT * retNode, //!< [out] To capture parent of returned bottom node.
     uint8_t & retIdx //!< [out] To capture sibling idx of returned bottom node.
     ) const noexcept {
      static_assert(ROW < ROW_NUM, "ROW should be smaller than ROW_NUM");
      assert(pos < this->getSumOfWeight<ROW>());

      retNode = this;
      while (true) {
        retIdx = 0;
        auto array = retNode->getConstPSumArray<ROW>();
        while (pos >= array[retIdx + 1]) {
          ++retIdx;
        }
        pos -= array[retIdx];
        if (isBorder()) {
          return;
        }
        retNode = retNode->getChildPtr(retIdx);
      }
    }


  private:
    //// private modifier (intend to call them from member function of BTreeNode)
    void setLmBorderNode(BTreeNodeT * lmBorderNode) noexcept {
      lmBorderNode_ = lmBorderNode;
    }


    void updateLmBorderNode(BTreeNodeT * lmBorderNode) noexcept {
      auto * node = this;
      while (true) {
        node->setLmBtm(lmBorderNode);
        if (node->isRoot() || node->getIdxInSibling() > 0) {
          break;
        }
        node = node->getParent();
      }
    }


    void unroot() noexcept {
      flags_ &= ~isRootBit;
    }


    void setParentRef(BTreeNodeT * newParent, uint8_t newIdxInSibling) noexcept {
      this->parent_ = newParent;
      this->idxInSibling_ = newIdxInSibling;
    }


    void setChildPtr(BTreeNodeT * child, uint8_t idx) noexcept {
      assert(idx < numChildren_);
      children_[idx] = child;
    }


    void makeNewRoot(BTreeNodeT * fstHalf, BTreeNodeT * sndHalf) {
      auto newRoot = new BTreeNodeT(false, true, fstHalf->getLmBorderNode());
      auto * parent = fstHalf->getParent();
      if (parent != NULL) { // BTrees are stacked
        const auto idxInSib = fstHalf->getIdxInSibling();
        parent->setChildPtr(newRoot, idxInSib); // parent points to newRoot
        newRoot->setParentRef(parent, idxInSib); // newRoot points to parent
      }
      newRoot->pushbackBTreeNode(fstHalf);
      newRoot->pushbackBTreeNode(sndHalf);
    }


    /*!
     * @brief Put node of ::BTreeNodeT type at idx.
     * @note
     *   - "psum_" is set using immediately preceeding psum value and weights of "child".
     *   - Links between "this" and "child" is fixed.
     *   - "numChildren_" is not incremented.
     */
    void putBTreeNode
    (
     BTreeNodeT * child,
     const uint8_t idx
     ) noexcept {
      assert(idx < kB);

      children_[idx] = child;
      for (uint8_t row = 0; row < ROW_NUM; ++row) {
        psum_[row][idx+1] = psum_[row][idx] + child->getSumOfWeight<row>();
      }
      child->setParentRef(this, idx);
    }


    /*!
     * @brief Pushback node of ::BTreeNodeT type.
     * @note
     *   - Links between "this" and "child" is fixed.
     *   - "numChildren_" is incremented.
     */
    void pushbackBTreeNode
    (
     BTreeNodeT * child
     ) noexcept {
      assert(numChildren_ < kB);

      putBTreeNode(child, numChildren_);
      ++numChildren_;
    }


    /*!
     * @brief Put bottom node other than ::BTreeNodeT type.
     * @note
     *   Link from "child" to "this" should be set outside this function (if "child" maintains it).
     */
    void putBtm
    (
     BTreeNodeT * child,
     const uint8_t idx,
     const uint64_t (&weights)[ROW_NUM]
     ) noexcept {
      assert(isBorder());
      assert(idx < kB);

      children_[idx] = child;
      for (uint8_t row = 0; row < ROW_NUM; ++row) {
        psum_[row][idx+1] = psum_[row][idx] + weights[row];
      }
    }


    /*!
     * @brief Pushback bottom node other than ::BTreeNodeT type.
     * @note
     *   Link from "child" to "this" should be set outside this function (if "child" maintains it).
     */
    void pushbackBtm
    (
     BTreeNodeT * child,
     const uint64_t (&weights)[ROW_NUM]
     ) noexcept {
      assert(isBorder());
      assert(numChildren_ < kB);

      putBtm(child, numChildren_, weights);
      ++numChildren_;
    }


    void shiftR
    (
     const uint64_t (&add_weights)[ROW_NUM],
     const uint8_t idx, //!< Beginning idx to shift.
     const uint8_t shift
     ) {
      assert(idx < numChildren_);
      assert(numChildren_ + shift <= kB);

      for (uint8_t i = numChildren_; idx < i; --i) {
        auto child = children_[i-1];
        children_[i + shift - 1] = child;
        child->idxInSibling_ = i + shift - 1;
      }
      for (uint8_t row = 0; row < ROW_NUM; ++row) {
        for (uint8_t i = numChildren_; idx < i; --i) {
          psum_[row][i + shift] = psum_[row][i] + add_weights[row];
        }
      }
      numChildren_ += shift;
    }


    void overflowToL
    (
     BTreeNodeT * lnode,
     BTreeNodeT * sndHalf,
     const uint8_t idx
     ) {
      assert(idx < numChildren_);

      const auto lnum = lnode->getNumChildren();
      const auto numToLeft = lnum/2 + kB/2 - lnum;
      const bool idxIsInLeft = idx < numToLeft;
      {
        const auto stop = (idxIsInLeft)? idx+1 : numToLeft;
        for (uint8_t i = 0; i < stop; ++i) {
          lnode->pushbackBTreeNode(children_[i]);
        }
      }
      if (idxIsInLeft) {
        lnode->pushbackBTreeNode(sndHalf);
        for (uint8_t i = idx + 1; i < numToLeft; ++i) {
          lnode->pushbackBTreeNode(children_[i]);
        }
      }

      {
        this->setLmBorderNode(children_[numToLeft]->getLmBorderNode());
        const auto stop = (idxIsInLeft)? numChildren_ : idx+1;
        for (uint8_t i = numToLeft; i < stop; ++i) {
          putBTreeNode(children_[i], i - numToLeft);
        }
      }
      if (!idxIsInLeft) {
        putBTreeNode(sndHalf, idx + 1 - numToLeft);
        for (uint8_t i = idx + 1; i < numChildren_; ++i) {
          putBTreeNode(children_[i], i - numToLeft + 1);
        }
      }
      numChildren_ -= (numToLeft - idxIsInLeft);
    }


    void overflowToR
    (
     BTreeNodeT * rnode,
     BTreeNodeT * sndHalf,
     const uint8_t idx
     ) {
      assert(idx < numChildren_);

      uint64_t temp_weights[ROW_NUM];

      const auto rnum = rnode->getNumChildren();
      const auto leftEnd = rnum/2 + kB/2;
      const bool idxIsInLeft = idx < leftEnd;
      auto shift = numChildren_ - leftEnd + !idxIsInLeft;
      if (rnum) {
        for (uint8_t row = 0; row < ROW_NUM; ++row) {
          temp_weights[row] = psum_[row][numChildren_] - psum_[row][leftEnd];
        }
        rnode->shiftR(temp_weights, 0, shift); // Shift elements in rnode
      }

      if (rnode->isBorder()) {
        rnode->setLmBorderNode(rnode);
      } else {
        rnode->setLmBorderNode(children_[leftEnd]->getLmBorderNode());
      }
      const auto stop = (idxIsInLeft)? numChildren_ : idx+1;
      for (uint8_t i = leftEnd; i < stop; ++i) {
        rnode->putBTreeNode(children_[i], i - leftEnd);
      }
      if (!idxIsInLeft) {
        rnode->putBTreeNode(sndHalf, idx + 1 - leftEnd);
        for (uint8_t i = idx + 1; i < numChildren_; ++i) {
          rnode->putBTreeNode(children_[i], i + 1 - leftEnd);
        }
      }

      numChildren_ = leftEnd;
      if (idxIsInLeft) {
        for (uint8_t row = 0; row < ROW_NUM; ++row) {
          psum_[row][idx] -= sndHalf->getSumOfWeight<row>();
          temp_weights[row] = 0;
        }
        shiftR(temp_weights, idx + 1, 1);
        putBTreeNode(sndHalf, idx + 1);
        ++numChildren_;
      }
    }


    BTreeNodeT * shiftR_Btm
    (
     const uint64_t (&add_weights)[ROW_NUM],
     const uint8_t idx, //!< Beginning idx to shift.
     const uint8_t shift
     ) {
      assert(idx < numChildren_);
      assert(numChildren_ + shift <= kB);

      for (uint8_t i = numChildren_; idx < i; --i) {
        children_[i + shift - 1] = children_[i-1];
      }
      for (uint8_t row = 0; row < ROW_NUM; ++row) {
        for (uint8_t i = numChildren_; idx < i; --i) {
          psum_[row][i + shift] = psum_[row][i] + add_weights[row];
        }
      }
      numChildren_ += shift;
    }


    void overflowToL_Btm
    (
     BTreeNodeT * lnode,
     BTreeNodeT * sndHalf,
     const uint64_t (&weights)[ROW_NUM],
     const uint8_t idx
     ) {
      assert(idx < numChildren_);

      uint64_t temp_weights[ROW_NUM];

      const auto lnum = lnode->getNumChildren();
      const auto numToLeft = lnum/2 + kB/2 - lnum;
      const bool idxIsInLeft = idx < numToLeft;
      {
        const auto stop = (idxIsInLeft)? idx : numToLeft;
        for (uint8_t i = 0; i < stop; ++i) {
          getWeightOfChild(i, temp_weights);
          lnode->pushbackBtm(children_[i], temp_weights);
        }
      }
      if (idxIsInLeft) {
        getWeightOfChild(idx, temp_weights);
        for (uint8_t row = 0; row < ROW_NUM; ++row) {
          temp_weights[row] -= weights[row];
        }
        lnode->pushbackBtm(children_[idx], temp_weights);
        lnode->pushbackBtm(sndHalf, weights);
        for (uint8_t i = idx + 1; i < numToLeft; ++i) {
          getWeightOfChild(i, temp_weights);
          lnode->pushbackBtm(children_[i], temp_weights);
        }
      }
      //
      {
        const auto stop = (idxIsInLeft)? numChildren_ : idx;
        for (uint8_t i = numToLeft; i < stop; ++i) {
          getWeightOfChild(i, temp_weights);
          putBtm(children_[i], i - numToLeft, temp_weights);
        }
      }
      if (!idxIsInLeft) {
        getWeightOfChild(idx, temp_weights);
        for (uint8_t row = 0; row < ROW_NUM; ++row) {
          temp_weights[row] -= weights[row];
        }
        putBtm(children_[idx], idx - numToLeft, temp_weights);
        putBtm(sndHalf, idx + 1 - numToLeft, weights);
        for (uint8_t i = idx + 1; i < numChildren_; ++i) {
          getWeightOfChild(i, temp_weights);
          putBtm(children_[i], i - numToLeft + 1, temp_weights);
        }
      }
      numChildren_ -= (numToLeft - idxIsInLeft);
    }


    void overflowToR_Btm
    (
     BTreeNodeT * rnode,
     BTreeNodeT * sndHalf,
     const uint64_t (&weights)[ROW_NUM],
     const uint8_t idx
     ) {
      assert(idx < numChildren_);

      uint64_t temp_weights[ROW_NUM];

      const auto rnum = rnode->getNumChildren();
      const auto leftEnd = rnum/2 + kB/2;
      const bool idxIsInLeft = idx < leftEnd;
      auto shift = numChildren_ - leftEnd + !idxIsInLeft;
      if (rnum) {
        for (uint8_t row = 0; row < ROW_NUM; ++row) {
          temp_weights[row] = psum_[row][numChildren_] - psum_[row][leftEnd];
        }
        rnode->shiftR_Btm(temp_weights, 0, shift); // Shift elements in rnode
      }

      {
        const auto stop = (idxIsInLeft)? numChildren_ : idx;
        for (uint8_t i = leftEnd; i < stop; ++i) {
          getWeightOfChild(i, temp_weights);
          rnode->putBtm(children_[i], i - leftEnd, temp_weights);
        }
      }
      if (!idxIsInLeft) {
        getWeightOfChild(idx, temp_weights);
        for (uint8_t row = 0; row < ROW_NUM; ++row) {
          temp_weights[row] -= weights[row];
        }
        rnode->putBtm(children_[idx], idx - leftEnd, temp_weights);
        rnode->putBtm(sndHalf, idx + 1 - leftEnd, weights);
        for (uint8_t i = idx + 1; i < numChildren_; ++i) {
          getWeightOfChild(i, temp_weights);
          rnode->putBtm(children_[i], i + 1 - leftEnd, temp_weights);
        }
      }
      //
      numChildren_ = leftEnd;
      if (idxIsInLeft) {
        for (uint8_t row = 0; row < ROW_NUM; ++row) {
          psum_[row][idx] -= weights[row];
          temp_weights[row] = 0;
        }
        shiftR_Btm(temp_weights, idx + 1, 1);
        putBtm(sndHalf, idx + 1, weights);
        ++numChildren_;
      }
    }


  public:
    //// public modifier
    /*!
     * @brief Handle the situation where 'children_[idx]' is split to 'children_[idx]' and 'sndHalf' when child node is ::BTreeNode type.
     */
    void handleSplitOfChild
    (
     BTreeNodeT * sndHalf,
     const uint8_t idx
     ) {
      assert(idx <= numChildren_);

      const auto end = numChildren_;
      if (end < kB) { // Easy case: Current node is not full.
        numChildren_ = idx;
        this->pushbackBTreeNode(children_[idx]);
        auto * pushC = sndHalf;
        for (uint8_t i = idx+1; i <= end; ++i) {
          auto * tmp = children_[i];
          this->pushbackBTreeNode(pushC);
          pushC = tmp;
        }
        return;
      }

      if (!isRoot()) {
        // Check siblings if they have space.
        if (idxInSibling_) { // Check previous sibling.
          auto lnode = parent_->getChildPtr(idxInSibling_ - 1);
          if (lnode->getNumChildren() < kB - 1) { // Previous sibling is not full. (-1 for easier implementation)
            overflowToL(lnode, sndHalf, idx);
            return lnode;
          }
        }
        if (idxInSibling_ + 1 < parent_->getNumChildren()) { // Check next sibling.
          auto rnode = parent_->getChildPtr(idxInSibling_ + 1);
          if (rnode->getNumChildren() < kB - 1) { // Next sibling is not full. (-1 for easier implementation)
            overflowToR(rnode, sndHalf, idx);
            return rnode;
          }
        }
      }

      { // this node has to be split
        auto newNode = new BTreeNodeT(this->isBorder(), false, nullptr);
        overflowToR(newNode, sndHalf, idx);
        if (!isRoot()) {
          parent_->handleSplitOfChild(newNode, idxInSibling_);
        } else {
          this->unroot();
          makeNewRoot(this, newNode);
        }
      }
    }


    /*!
     * @brief Handle the situation where 'children_[idx]' is split to 'children_[idx]' and 'sndHalf' when child node is not ::BTreeNode type.
     */
    BTreeNodeT * handleSplitOfBtm
    (
     BTreeNodeT * sndHalf,
     const uint64_t (&weights)[ROW_NUM],
     const uint8_t idx
     ) {
      assert(isBorder());
      assert(idx <= numChildren_);

      if (numChildren_ < kB) { // Easy case: Current node is not full.
        uint64_t temp_weights[ROW_NUM] = {}; // 0 initialized.
        for (uint8_t row = 0; row < ROW_NUM; ++row) {
          psum_[row][idx] -= weights[row];
        }
        shiftR_Btm(temp_weights, idx + 1, 1);
        putBtm(sndHalf, idx + 1, weights);
        return NULL;
      }

      if (!isRoot()) {
        // Check siblings if they have space.
        if (idxInSibling_) { // Check previous sibling.
          auto lnode = parent_->getChildPtr(idxInSibling_ - 1);
          if (lnode->getNumChildren() < kB - 1) { // Previous sibling is not full. (-1 for easier implementation)
            overflowToL_Btm(lnode, sndHalf, weights, idx);
            return lnode;
          }
        }
        if (idxInSibling_ + 1 < parent_->getNumChildren()) { // Check next sibling.
          auto rnode = parent_->getChildPtr(idxInSibling_ + 1);
          if (rnode->getNumChildren() < kB - 1) { // Next sibling is not full. (-1 for easier implementation)
            overflowToR_Btm(rnode, sndHalf, weights, idx);
            return rnode;
          }
        }
      }

      { // This node has to be split.
        auto newNode = new BTreeNode(true, false, children_[kB/2]);
        overflowToR_Btm(newNode, sndHalf, weights, idx);
        if (!isRoot()) {
          parent_->handleSplitOfChild(newNode, idxInSibling_);
        } else {
          this->unroot();
          makeNewRoot(this, newNode);
        }
        return newNode;
      }
    }


    /*!
     * @brief Change weight of 'idx'-th child of this node, and accordingly all 'psum_' values needed to be fixed.
     */
    template<uint8_t row>
    void changePSumFrom
    (
     const uint8_t idx,
     const int64_t change
     ) noexcept {
      assert(row < ROW_NUM);

      for (uint8_t i = idx; i < numChildren_; ++i) {
        psum_[row][i] += change;
      }
      if (parent_ != NULL) { // we do not use isRoot() here for convenience. That is, when we stack two or more BTrees, the change will be propagated.
        parent_->changePSumFrom(row, idxInSibling_, change);
      }
    }


  public:
    //// calculate statistics
    size_t calcMemBytes() const noexcept {
      size_t sumOfSize = sizeof(*this);
      if (!isBorder()) {
        for (uint8_t i = 0; i < numChildren_; ++i) {
          sumOfSize += children_[i]->calcMemBytes();
        }
      }
      return sumOfSize;
    }


    size_t calcNumUsed() const noexcept {
      size_t numOfUsed = numChildren_;
      if (!isBorder()) {
        for (uint8_t i = 0; i < numChildren_; ++i) {
          numOfUsed += children_[i]->calcNumUsed();
        }
      }
      return numOfUsed;
    }


    size_t calcNumSlots() const noexcept {
      size_t numOfSlots = kB;
      if (!isBorder()) {
        for (uint8_t i = 0; i < numChildren_; ++i) {
          numOfSlots += children_[i]->calcNumSlots();
        }
      }
      return numOfSlots;
    }
  };
} // namespace itmmti

#endif
