/*!
 * Copyright (c) 2017 Tomohiro I
 *
 * This program is released under the MIT License.
 * http://opensource.org/licenses/mit-license.php
 */
/*!
 * @file PSumWithValue.hpp
 * @brief Dynamic partial sum data structure implemented by B+tree, where each leaf is associated with value.
 * @author Tomohiro I
 * @date 2017-10-27
 */
#ifndef INCLUDE_GUARD_BTreeWithValue
#define INCLUDE_GUARD_BTreeWithValue

#include <stdint.h>

#include <iostream>
#include <algorithm>
#include <cassert>
#include <string>
#include <fstream>
#include <sstream>

// include from Basics
#include "BitsUtil.hpp"
#include "StepCode.hpp"
#include "MemUtil.hpp"
//

#include "BTree.hpp"


namespace itmmti
{
  template <uint8_t kB, uint8_t ROW_NUM> class BTreeNode;


  template <uint16_t kBtmB>
  class BtmNodeForPSumWithVal
  {
  private:
    //// Private constant, alias etc.
    using BtmNodeT = BtmNodeForPSumWithVal<kBtmB>;


  private:
    //// Private member variables.
    StepCodeCore<2 * kBtmB> stcc_; //!< weights and values are stored alternatingly.
    uint16_t bitCapacity_; //!< Current bit capacity.
    uint16_t bitSize_; //!< Current bit size.
    uint8_t size_; //!< Current size (number of elements).
    uint8_t wCodesAuxM_[2 * kBtmB / StepCodeUtil::kWCNum - 1];


  public:
    BtmNodeForPSumWithVal
    (
     uint16_t initBitCapa = 0;
     ) : stcc_(), bitCapacity_(0), bitSize_(0), size_(0)
    {
      stcc_.changeBitCapacity(bitCapacity_, 0, initBitCapa);
    }


    ~BtmNodeForPSumWithVal()
    {
      // stcc_ is freed.
    }


    uint16_t getNumChildren() const noexcept
    {
      return size_;
    }


    uint64_t getWeight
    (
     uint16_t childIdx
     ) const noexcept {
      assert(childIdx < static_cast<uint16_t>(size_));

      return stcc_.read(2 * childIdx);
    }


    uint64_t getVal
    (
     uint16_t childIdx
     ) const noexcept {
      assert(childIdx < static_cast<uint16_t>(size_));

      return stcc_.read(2 * childIdx + 1);
    }


    uint64_t getSumOfWeight() const noexcept
    {
      uint64_t sum = 0;
      uint64_t bitPos = 0;
      for (uint64_t i = 0; i < size_; ++i) {
        const auto w = stcc_.readW(2 * i);
        sum += stcc_.readWBits(bitPos, w, bits::UINTW_MAX(w));
        bitPos += w + stcc_.readW(2 * i + 1);
      }

      return sum;
    }


    /*!
     * @brief Calculate the beginning bit-pos of "idx"-th value in stcc_.
     */
    uint64_t calcBitPos
    (
     const uint16_t idx, //!< in [0, 2*size_]
     ) const noexcept {
      assert(idx <= 2 * static_cast<uint16_t>(size_));

      return stcc_.calcBitPos(idx, wCodesAuxM_);
    }


    /*!
     * @brief Return child index where partial sum of "pos+1" is achieved.
     * @pre Answer should be found in this bottom node.
     */
    uint16_t searchPos
    (
     uint64_t & pos, //!< [in,out] Give position to search. It is modified to relative position.
     uint64_t & retWeight, //!< [out] Capture last weight where pos exists.
     uint64_t & retVal //!< [out] Capture value.
     ) const noexcept {
      assert(pos < getSumOfWeight()); // pre

      uint16_t bitPos = 0;
      uint16_t i = 0;
      while (true) {
        const auto w = stcc_.readW(i);
        retWeight = stcc_.readWBits(bitPos, w, bits::UINTW_MAX(w));
        const auto w_val = stcc_.readW(++i);
        if (pos < retWeight) {
          retVal = stcc_.readWBits(bitPos + w, w_val, bits::UINTW_MAX(w_val));
          return i / 2;
        }
        bitPos += w + w_val;
        pos -= retWeight;
        ++i;
      }
    }


    /*!
     * @brief Return child index where partial sum of "pos+1" is achieved.
     * @pre Answer should be found in this bottom node.
     */
    uint16_t searchPos
    (
     uint64_t & pos, //!< [in,out] Give position to search. It is modified to relative position.
     uint64_t & retWeight //!< [out] Capture last weight where pos exists.
     ) const noexcept {
      assert(pos < getSumOfWeight()); // pre

      uint16_t bitPos = 0;
      uint16_t i = 0;
      while (true) {
        const auto w = stcc_.readW(i);
        retWeight = stcc_.readWBits(bitPos, w, bits::UINTW_MAX(w));
        if (pos < retWeight) {
          return i/2;
        }
        bitPos += w + stcc_.readW(++i);
        pos -= retWeight;
        ++i;
      }
    }


    /*!
     * @brief Return child index where partial sum of "pos+1" is achieved.
     * @pre Answer should be found in this bottom node.
     */
    uint16_t searchPos
    (
     uint64_t & pos //!< [in,out] Give position to search. It is modified to relative position.
     ) const noexcept {
      assert(pos < getSumOfWeight()); // pre

      uint64_t weight;
      return searchPos(pos, weight);
    }


  private:
    /*!
     * @brief Get read-only array pointer.
     */
    const uint64_t * getConstPtr_vals() const noexcept
    {
      return stcc_.getConstPtr_vals();
    }


    /*!
     * @brief Get read-only wCodes_ array pointer.
     */
    const uint64_t * getConstPtr_wCodes() const noexcept
    {
      return stcc_.getConstPtr_wCodes();
    }


    /*!
     * @brief Resize "size_" to "newSize".
     * @note
     *   It does not change bitCapacity.
     */
    void resize
    (
     const uint8_t newSize
     ) noexcept {
      assert(newSize <= kBtmB);

      size_ = newSize;
    }


    void changeWCodes
    (
     const uint64_t * srcWCodes,
     const uint16_t srcIdxBeg, //!< Beginning idx of src.
     const uint16_t srcLen, //!< Length of wCodes of src to insert.
     const uint16_t tgtIdxBeg, //!< Beginning idx of tgt.
     const uint16_t tgtLen, //!< Length of wCodes of tgt to delete.
     ) noexcept {
      assert(tgtIdxBeg + tgtLen <= static_cast<uint16_t>(size_) * 2);
      assert(tgtIdxBeg - tgtLen + srcLen <= kBtmB * 2);

      uint16_t sizeTimes2 = static_cast<uint16_t>(size_) * 2;
      stcc_.changeWCodes(srcWCodes, srcIdxBeg, srcLen, tgtIdxBeg, tgtLen, sizeTimes2);
      this->resize(sizeTimes2 / 2);
    }


    /*!
     * @brief Change (shift) value positions.
     * @pre The resulting size should be within capacity.
     */
    void changeValPos
    (
     const uint64_t bitPos, //!< bitPos of bv_.
     const uint64_t insBitLen, //!< Bit-length to insert in bv_
     const uint64_t delBitLen //!< Bit-length to delete in bv_
     ) noexcept {
      stcc_.changeValPos(bitSize_, bitPos, insBitLen, delBitLen);
    }


    /*!
     * @brief Move vals.
     * @pre Bit-regions should be within capacity.
     */
    void mvVals
    (
     const uint64_t * srcVals,
     const uint64_t srcBitPos,
     const uint64_t tgtBitPos,
     const uint64_t bitLen
     ) noexcept {
      stcc_.mvVals(srcVals, srcBitPos, tgtBitPos, bitLen);
    }


    /*!
     * @brief update wCodesAuxM.
     */
    void updateWCodesAuxM
    (
     const uint16_t idxBeg,
     const uint16_t idxEnd,
     ) noexcept {
      assert(idxBeg < idxEnd);

      const uint64_t beg = (idxBeg - 1) / StepCodeUtil::kWCNum;
      const uint64_t end = (idxEnd - 1) / StepCodeUtil::kWCNum + (idxEnd <= (kMaxNum - StepCodeUtil::kWCNum));
      stcc_.updateWCodesAuxM(wCodesAuxM_, beg, end);
    }


    void overflowToL_wCodes
    (
     BtmNodeT * lnode,
     const uint64_t * srcWCodes,
     const uint16_t numChild_ins,
     const uint16_t childIdx, //!< Beginning childIdx of tgt.
     const uint16_t numChild_del //!< Length of wCodes of tgt to delete.
     ) noexcept {
      assert(childIdx < numChildren_);

      std::tuple<uint64_t, uint64_t, uint64_t> changeList[2];
      uint8_t clSize = 0;
      const uint16_t numL_old = lnode->getNumChildren();
      const uint16_t numR_old = this->getNumChildren();
      const uint16_t numTotal = numL_old + numR_old + numChild_ins - numChild_del;
      const uint16_t numL_new = num_total/2;
      const uint16_t numR_new = num_total - numL_new;
      const uint16_t numToLeft = numL_new - numL_old;
      const bool isNewElemInL = childIdx < numToLeft;
      const bool isNewElemInR = childIdx + numChild_ins > numToLeft;
      uint64_t sumWL = lnode->bitSize_;
      uint16_t sizeTimes2 = numL_old * 2;
      uint16_t curNumAfterDel = 0;
      uint16_t curNumSrcWCodes = 0;
      {
        const auto num = (isNewElemInL)? childIdx : numToLeft;
        if (num) {
          const auto w = this->calcBitPos(2 * num);
          changeList[clSize++] = {0, sumWL, w};
          sumWL += w;
          lnode_->changeWCodes(this->getConstPtr_wCodes(), 0, 2 * num, sizeTimes2, 0, sizeTimes2);
        }
      }
      if (isNewElemInL) {
        curNumSrcWCodes = std::min(numToLeft - childIdx, numChild_ins);
        sumWL += StepCodeUtil::sumW(srcWCodes, 0, 2 * curNumSrcWCodes);
        lnode_->stcc_.changeWCodes(srcWCodes, 0, 2 * curNumSrcWCodes, sizeTimes2, 0, sizeTimes2);
        if (numToLeft > childIdx + numChild_ins) {
          curNumAfterDel = numToLeft - (childIdx + numChild_ins);
          const auto w = this->stcc_.sumW(2 * (childIdx + numChild_del), 2 * (childIdx + numChild_del + curNumAfterDel));
          changeList[clSize++] = {this->calcBitPos(2 * (childIdx + numChild_del)), sumWL, w};
          sumWL += w;
          lnode_->stcc_.changeWCodes(this->getConstPtr_wCodes(), 2 * (childIdx + numChild_del), 2 * curNumAfterDel, sizeTimes2, 0, sizeTimes2);
        }
      }
      lnode->size_ = numL_new;
      lnode->updateWCodesAuxM(lnode->wCodesAuxM_, 2 * numL_old, 2 * numL_new);
      { // Update vals of lnode.
        if (lnode->bitCapacity_ < sumWL) {
          lnode->stcc_.changeBitCapacity(lnode->bitCapacity_, lnode->bitSize_, sumWL);
        }
        for (uint8_t i = 0; i < clSize; ++i) {
          lnode->mvVals(this->getConstPtr_vals(), std::get<0>(changeList[i]), std::get<1>(changeList[i]), std::get<2>(changeList[i]));
        }
        lnode->bitSize_ = sumWL;
      }

      // Update this node.
      clSize = 0;
      sizeTimes2 = numR_old * 2;
      const auto bitPosOfLastChunk = this->calcBitPos(2 * (childIdx + numChild_del + curNumAfterDel));
      const auto bitSizeOfLastChunk = this->bitSize_ - bitPosOfLastChunk;
      uint64_t sumWR = bitSizeOfLastChunk;
      if (numToLeft < childIdx) {
        const uint16_t num = childIdx - numToLeft;
        const uint16_t bitPos_src = this->calcBitPos(2 * numToLeft);
        const uint8_t w = this->calcBitPos(2 * childIdx) - bitPos_src;
        sumWR += w;
        changeList[clSize++] = {bitPos_src, 0, w};
        this->stcc_.changeWCodes(this->getConstPtr_wCodes(), 2 * numToLeft, 2 * num, 0, 2 * num, sizeTimes2);
      }
      if (isNewElemInR) {
        sumWR += StepCodeUtil::sumW(srcWCodes, 2 * curNumSrcWCodes, 2 * numChild_ins);
      }
      if (numR_old != numR_new && numR_old != childIdx + numChild_del) { // Need shift wCodes of "this" node.
        const uint16_t srcBeg = childIdx + numChild_del + curNumAfterDel;
        const uint16_t num = numR_old - srcBeg;
        const uint16_t tgtBeg = numR_new - num;
        this->changeWCodes(this->getConstPtr_wCodes(), 2 * srcBeg, 2 * num, 2 * tgtBeg, 2 * num, sizeTimes2);
        changeList[clSize++] = {bitPosOfLastChunk, sumWR - bitSizeOfLastChunk, bitSizeOfLastChunk};
      }
      if (isNewElemInR) {
        const uint16_t num = numChild_ins - curNumSrcWCodes;
        this->changeWCodes(srcWCodes, 2*curNumSrcWCodes, 2*num, 2*(childIdx + curNumSrcWCodes - numToLeft), 2*num, sizeTimes2);
      }
      this->size_ = numR_new;
      this->updateWCodesAuxM(0, 2 * numR_new);
      { // Update vals of "this" node.
        if (this->bitCapacity_ < sumWR) {
          this->stcc_.changeBitCapacity(this->bitCapacity_, this->bitSize_, sumWR);
        }
        for (uint8_t i = 0; i < clSize; ++i) {
          this->mvVals(this->getConstPtr_vals(), std::get<0>(changeList[i]), std::get<1>(changeList[i]), std::get<2>(changeList[i]));
        }
        this->bitSize_ = sumWR;
      }
    }


    void overflowToR_wCodes
    (
     BtmNodeT * rnode,
     const uint64_t * srcWCodes,
     const uint16_t numChild_ins,
     const uint16_t childIdx, //!< Beginning childIdx of tgt.
     const uint16_t numChild_del //!< Length of wCodes of tgt to delete.
     ) noexcept {
      assert(childIdx < numChildren_);

      std::tuple<uint64_t, uint64_t, uint64_t> changeList[2];
      uint8_t clSize = 0;
      const uint16_t numL_old = this->getNumChildren();
      const uint16_t numR_old = rnode->getNumChildren();
      const uint16_t numTotal = numL_old + numR_old + numChild_ins - numChild_del;
      const uint16_t numL_new = num_total/2;
      const uint16_t numR_new = num_total - numL_new;
      const uint16_t numToRight = numR_new - numR_old;
      uint16_t sizeTimes2 = 2 * kBtmB; // dummy.

      uint16_t numSrcWCodesInL = 0;
      uint16_t numToRight1 = 0;
      uint16_t numToRight2 = 0;
      if (childIdx < numL_new) { // new elements are in L
        if (childIdx + numChild_ins <= numL_new) { // new elements are only in L
          numSrcWCodesInL = numChild_ins;
          numToRight2 = numToRight;
        } else { // new elements are also in in R
          numSrcWCodesInL = numL_new - childIdx;
          numToRight2 = numR_old - (childIdx + numChild_del);
        }
      } else { // new elements are in R
        numToRight1 = childIdx - numL_new;
        numToRight2 = numL_old - (childIdx + numChild_del);
      }

      if (numR_old) { // shift wCodes of R to make space
        rnode->stcc_.changeWCodes(rnode->getConstPtr_wCodes(), 0, 2 * numToRight, 2 * numToRight, 2 * numToRight, sizeTimes2);
      }

      uint16_t numR_inc = 0;
      uint64_t sumWR_inc = 0;
      if (numToRight1) {
        const auto bitPos = this->calcBitPos(2 * (childIdx - numToRight1));
        const auto w = this->calcBitPos(2 * childIdx) - bitPos;
        changeList[clSize++] = {bitPos, 0, w};
        sumWR_inc += w;
        rnode->changeWCodes(this->getConstPtr_wCodes(), 2 * (childIdx - numToRight1), 2 * numToRight1, 0, 2 * numToRight1, sizeTimes2);
        numR_inc += numToRight1;
      }
      if (numSrcWCodesInL != numChild_ins) {
        sumWR_inc += StepCodeUtil::sumW(srcWCodes, 2 * numSrcWCodesInL, 2 * numChild_ins);
        rnode->changeWCodes(srcWCodes, 2 * numSrcWCodesInL, 2 * (numChild_ins - numSrcWCodesInL), numR_inc, 2 * (numChild_ins - numSrcWCodesInL), sizeTimes2);
        numR_inc += numSrcWCodesInL;
      }
      if (numToRight2) {
        const auto bitPos = this->calcBitPos(2 * (numL_old - numToRight2));
        const auto w = this->bitSize_ - bitPos;
        changeList[clSize++] = {bitPos, sumWR_inc, w};
        sumWR_inc += w;
        rnode->changeWCodes(this->getConstPtr_wCodes(), 2 * (numL_old - numToRight2), 2 * numToRight2, numR_inc, 2 * numToRight2, sizeTimes2);
      }
      rnode->size_ = numR_new;
      rnode->updateWCodesAuxM(0, 2 * numR_new);
      { // Update vals of "rnode".
        if (rnode->bitCapacity_ < rnode->bitSize_ + sumWR_inc) {
          rnode->stcc_.changeBitCapacity(rnode->bitCapacity_, rnode->bitSize_, rnode->bitSize_ + sumWR_inc);
        }
        if (numR_old) {
          rnode->mvVals(rnode->getConstPtr_vals(), 0, sumWR_inc, rnode->bitSize_);
        }
        for (uint8_t i = 0; i < clSize; ++i) {
          rnode->mvVals(this->getConstPtr_vals(), std::get<0>(changeList[i]), std::get<1>(changeList[i]), std::get<2>(changeList[i]));
        }
        rnode->bitSize_ += sumWR_inc;
      }

      if (numSrcWCodesInL) {
        const auto sumW_insert = StepCodeUtil::sumW(srcWCodes, 0, 2 * numSrcWCodesInL);
        const auto numTail = numL_new - (childIdx + numChild_ins);
        if (numTail) {
          const auto bitPos = this->calcBitPos(2 * (childIdx + numChild_del));
          const auto w = this->calcBitPos(2 * (childIdx + numChild_del + numTail)) - bitPos;
          const auto newBitSize = bitPos + sumW_insert + w;
          if (this->bitCapacity_ < newBitSize) {
            this->stcc_.changeBitCapacity(this->bitCapacity_, this->bitSize_, newBitSize);
          }
          this->mvVals(this->getConstPtr_vals(), bitPos, bitPos + sumW_insert, w);
          this->bitSize_ = newBitSize;
        }
        sizeTimes2 = 2 * (childIdx + numChild_del + numTail);
        this->changeWCodes(srcWCodes, 0, 2 * numSrcWCodesInL, 2 * childIdx, 2 * numChild_del, sizeTimes2);
        this->updateWCodesAuxM(2 * childIdx, 2 * numR_new);
      } else {
        this->bitSize_ = this->calcBitPos(2 * numL_new); // shrink (just change bitSize)
      }
      this->size_ = numL_new;
    }


  public:
    void insert
    (
     const uint64_t * weights,
     const uint64_t * vals,
     const uint16_t numChild_ins,
     const uint16_t childIdx, //!< Beginning childIdx of tgt.
     const uint16_t numChild_del, //!< Length of wCodes of tgt to delete.
     BTreeNodeT * parent,
     const uint8_t idxInSibling
     ) noexcept {
      assert(numChild_ins <= size);
      assert(childIdx + numChild_del <= size);

      uint64_t wCodesTemp[2 * kBtmB / StepCodeUtil::kWCNum];
      uint16_t sumW_ins = 0;
      uint64_t sumWeights_ins = 0;
      for (uint16_t ci = 0; ci < numChild_ins; ++i) {
        sumWeights_ins = weights[ci];
        uint8_t w = StepCodeUtil::calcSteppedW(weights[ci]);
        sumW_ins += w;
        StepCodeUtil::writeWCode(StepCodeUtil::calcWCodeFromSteppedW(w), wCodesTemp, 2 * ci);
        w = StepCodeUtil::calcSteppedW(vals[ci]);
        sumW_ins += w;
        StepCodeUtil::writeWCode(StepCodeUtil::calcWCodeFromSteppedW(w), wCodesTemp, 2 * ci + 1);
      }

      uint16_t sumW_del = 0;
      uint64_t sumWeights_del = 0;
      if (numChild_del) {
        uint64_t bitPos = this->calcBitPos(2*childIdx);
        sumW_del = bitPos; // Temporarily (later fixed by "sumW_del = bitPos - sumW_del").
        for (uint16_t ci = childIdx; ci < childIdx + numChild_del; ++ci) {
          const uint8_t w = stcc_.readW(2 * ci);
          sumWeights_del += stcc_.read(bitPos, w, bits::UINTW_MAX(w));
          bitPos += w + stcc_.readW(2 * ci + 1);
        }
        sumW_del = bitPos - sumW_del;
      }

      if (sumWeights_ins != sumWeights_del) {
        parent->changePSumFrom<ROW0>(idxInSibling, static_cast<int64_t>(sumWeights_ins) - sumWeights_del);
      }

      BtmNodeT * lnode = this;
      BtmNodeT * rnode = this;
      uint16_t childIdx_ins = childIdx;

      if (childIdx - numChild_del + numChild_ins <= kBtmB) { // Easy case: This node can accommodate inserting elements.
        uint16_t sizeTimes2 = static_cast<uint16_t>(size_) * 2;
        stcc_.changeWCodes(wCodesTemp, 0, 2 * numChild_ins, 2 * childIdx, 2 * numChild_del, sizeTimes2);
        size_ = sizeTimes2 / 2;
        this->updateWCodesAuxM(2 * childIdx, sizeTimes2);
        if (bitCapacity_ < bitSize_ + sumW_ins - sumW_del) { // Need expand bits.
          stcc_.changeBitCapacity(bitCapacity_, bitSize_, bitSize_ + sumW_ins - sumW_del);
        }
        goto insertNewElem;
      }

      if (idxInSibling) { // Check previous sibling.
        lnode = reinterpret_cast<BtmNodeT *>(parent_->getChildPtr(idxInSibling_ - 1));
        if (lnode->getNumChildren() <= kBtmB + numChild_del - numChild_ins) { // Previous sibling can accommodate overflowed elements.
          childIdx_ins += lnode->getNumChildren();
          this->overflowToL_wCodes(lnode, wCodesTemp, numChild_ins, childIdx, numChild_del);
          goto insertNewElem;
        }
      }
      lnode = this;
      if (idxInSibling_ + 1 < parent_->getNumChildren()) { // Check next sibling.
        rnode = reinterpret_cast<BtmNodeT *>(parent_->getChildPtr(idxInSibling_ + 1));
        if (rnode->getNumChildren() <= kBtmB + numChild_del - numChild_ins) { // Next sibling can accommodate overflowed elements.
          this->overflowToR(rnode, wCodesTemp, numChild_ins, childIdx, numChild_del);
          goto insertNewElem;
        }
      }

      { // This node has to be split
        rnode = new BtmNodeT();
        this->overflowToR_wCodes(rnode, wCodesTemp, numChild_ins, childIdx, numChild_del);
        const uint64_t weightOfRight[] = {rnode->getSumOfWeight()};
        parent->handleSplitOfBtm(reinterpret_cast<BTreeNodeT *>(newNode), weightOfRight, idxInSibling);
      }

    insertNewElem:
      uint16_t numL = lnode->getNumChildren();
      uint16_t numNewElemInL = 0
        if (childIdx_ins < numL) {
          uint64_t bitPos = lnode->calcBitPos(2 * childIdx);
          numNewElemInL = std::min(numChild_ins, numL - childIdx_ins);
          for (uint16_t ci = childIdx_ins; ci < childIdx_ins + numNewElemInL; ++ci) {
            uint8_t w = stcc_.readW(2 * ci);
            lnode->stcc_.writeWBits(weights[ci - childIdx_ins], bitPos, w);
            bitPos += w;
            w = stcc_.readW(2 * ci + 1);
            lnode->stcc_.writeWBits(vals[ci - childIdx_ins], bitPos, w);
            bitPos += w;
          }
        }
      if (numNewElemInL < numChild_ins) {
        childIdx_ins += numNewElemInL - numL;
        for (uint16_t ci = childIdx_ins; ci < childIdx_ins + numChild_ins - numNewElemInL; ++ci) {
          uint8_t w = stcc_.readW(2 * ci);
          rnode->stcc_.writeWBits(weights[ci - childIdx_ins], bitPos, w);
          bitPos += w;
          w = stcc_.readW(2 * ci + 1);
          rnode->stcc_.writeWBits(vals[ci - childIdx_ins], bitPos, w);
          bitPos += w;
        }
      }
    }


    /*!
     * @brief Replace weights and/or values.
     */
    void replace
    (
     const uint64_t * array, //!< Storing weights and vals alternatingly (starting with weight iff "idx = 0 mod 2").
     const uint16_t num, //!< Number of elements to replace.
     const uint16_t idx, //!< in [0..2*size_]. Beginning idx of tgt.
     BTreeNodeT * parent,
     const uint8_t idxInSibling
     ) {
      assert(idx + num <= 2 * static_cast<uint16_t>(size_));

      const uint16_t bitPos0 = this->calcBitPos(idx);
      uint16_t bitPos = bitPos0;
      uint16_t sumW_ins = 0;
      uint16_t sumW_del = 0;
      int64_t weights_diff = 0;
      for (uint16_t i = idx; i < idx + num; ++i) {
        const uint8_t w_old = stcc_.readW(i);
        const uint8_t w_new = StepCodeUtil::calcSteppedW(array[i - idx]);
        sumW_del += w_old;
        sumW_ins += w_new;
        if (i / 2) {
          weights_diff += array[i - idx] - stcc_.readWBits(bitPos, w_old);
        }
        stcc_.writeWCode(StepCodeUtil::calcWCodeFromSteppedW(w_new), i);
        bitPos += w_old;
      }
      this->updateWCodesAuxM(idx, idx + num);

      if (sumW_ins != sumW_del) {
        if (sumW_ins > sumW_del && bitCapacity_ > bitSize_ + sumW_ins - sumW_del) {
          stcc_.changeBitCapacity(bitCapacity_, bitSize_, bitSize_ + sumW_ins - sumW_del);
        }
        this->changeValPos(bitPos0, sumW_ins, sumW_del);
      }
      bitPos = bitPos0;
      for (uint16_t i = idx; i < idx + num; ++i) {
        uint8_t w = stcc_.readW(2 * ci);
        stcc_.writeWBits(array[i - idx], bitPos, w);
        bitPos += w;
      }

      if (weights_diff != 0) {
        parent->changePSumFrom<kRow0>(idxInSibling, weights_diff);
      }
    }


  public:
    //// statistics
    size_t calcMemBytes
    (
     bool includeThis = true
     ) const noexcept {
      size_t size = sizeof(*this) * includeThis;
      return size + bitCapacity_ / 8;
    }


    size_t calcMemBytesDynArray() const noexcept {
      return bitCapacity_ / 8;
    }
  };




  /*!
   * @brief Dynamic partial sum data structure implemented by B+tree, where each leaf is associated with value.
   * @tparam B Arity for internal node of B+tree, which should be in {32, 64, 128}.
   * @tparam BL Arity for bottom node of B+tree, which should be in {32, 64, 128}.
   * @par Notation
   *   - T: Current string represented by RLE.
   *   - Mixed tree: B+tree representing RLE of T.
   *     - btmM: Index of bottom node of B+tree on the array-based implementation.
   *             Each bottom node 'btmM' can have 'B' children, which correspond to indexes [btmM * B, (btmM+1) * B).
   *     - idxM: Indexes that are corresponding to children of btmM.
   *   - Separated tree: B+tree separately representing runs for each character.
   *     - btmS: Index of bottom node of B+tree on the array-based implementation (all separated trees share arrays).
   *             Each bottom node 'btmS' can have 'B' children, which correspond to indexes [btmS * B, (btmS+1) * B).
   *     - idxS: Indexes that are corresponding to children of btmS.
   */
  template <uint8_t kB = 32, uint8_t kBtmB = 32>
  class PSumWithValue
  {
  private:
    // Private constant, alias etc.
    using BTreeNodeT = BTreeNode<kB, kBtmB>;
    using BtmNodeT = BtmNodeForPSumWithVal<kBtmB>;


  public:
    // Public constant, alias etc.
    static constexpr uint8_t kRow0{0};


  private:
    // Private member variables.
    BTreeNodeT * root_; //!< Root of B+tree.


  public:
    PSumWithValue() :
      root_(nullptr)
    {}


    ~PSumWithValue() {
      clearAll();
    }


    /*!
     * @brief Initialize data structure with empty elements.
     */
    void init() {
      if (isReady()) {
        clearAll();
      }
      BtmNodeT * fstBtmNode = new BtmNodeT(kBtmB * 8);
      root_ = new BTreeNodeT(true, true, reinterpret_cast<BTreeNodeT *>(fstBtmNode));
      root_->pushbackBtm(reinterpret_cast<BTreeNodeT *>(fstBtmNode), {0});
      // Insert sentinel.
      const uint64_t dummyArray[] = {0};
      fstBtmNode->insert(dummyArray, dummyArray, 1, 0, 0, root_, 0);
    }


    /*!
     * @brief Free/delete all allocated objects.
     */
    void clearAll() {
      if (!isReady()) { // already cleared
        return;
      }
      // Delete bottom nodes.
      uint8_t idxInSib = 0;
      auto borderNode = root_->getLmBorderNode();
      while (reinterpret_cast<uintptr_t>(borderNode) != BTreeNodeT::NOTFOUND) {
        delete reinterpret_cast<BtmNodeT *>(borderNode->getChildPtr(idxInSib));
        borderNode = getNextBtmRef(idxInSib);
      }
      // Delete uppdar nodes.
      root_->clearUpperNodes();
      root_ = nullptr;
    }


    /*!
     * @brief Return if data structure is ready.
     */
    bool isReady() const noexcept {
      return (root_ != NULL);
    }


    /*!
     * @brief Return |T|.
     */
    size_t getSumOfWeight() const noexcept {
      assert(isReady());

      return root_->getSumOfWeight<kRow0>();
    }


  public:
    //// Search functions
    /*!
     * @brief Search for bottom node achieving psum of "pos+1".
     * @attention "pos" is modified to be the relative position (0base) from beginning of bottom.
     */
    void searchBtm
    (
     uint64_t & pos //!< [in,out] Give position to search (< |T|). It is modified to relative position.
     BTreeNodeT * retNode, //!< [out] To capture parent of returned bottom node.
     uint8_t & retIdx //!< [out] To capture sibling idx of returned bottom node.
     ) const noexcept {
      assert(isReady());
      assert(pos < root_->getSumOfWeight<kRow0>());

      root_->searchPos(pos, retNode, retIdx);
    }


  public:
    //// statistics
    size_t calcMemBytesUpperPart() const noexcept {
      return root_->calcMemBytes();
    }


    size_t calcMemBytesBtmPart() const noexcept {
      size_t size = 0;
      uint8_t idxInSib = 0;
      auto borderNode = root_->getLmBorderNode();
      while (reinterpret_cast<uintptr_t>(borderNode) != BTreeNodeT::NOTFOUND) {
        size += reinterpret_cast<BtmNodeT *>(borderNode->getChildPtr(idxInSib))->calcMemBytes();
        borderNode = getNextBtmRef(idxInSib);
      }
      return size;
    }


    size_t calcMemBytesDynArrayOfStepCode() const noexcept {
      size_t size = 0;
      uint8_t idxInSib = 0;
      auto borderNode = root_->getLmBorderNode();
      while (reinterpret_cast<uintptr_t>(borderNode) != BTreeNodeT::NOTFOUND) {
        size += reinterpret_cast<BtmNodeT *>(borderNode->getChildPtr(idxInSib))->calcMemBytesDynArray();
        borderNode = getNextBtmRef(idxInSib);
      }
      return size;
    }


    size_t calcMemBytes
    (
     bool includeThis = true
     ) const noexcept {
      size_t size = sizeof(*this) * includeThis;
      size += calcMemBytesUpperPart();
      size += calcMemBytesBtmPart();
      return size;
    }


    size_t calcNumUsedUpperPart() const noexcept {
      return root_->calcNumUsed();
    }


    size_t calcNumSlotsUpperPart() const noexcept {
      return root_->calcNumSlots();
    }


    size_t calcNumUsedBtmPart() const noexcept {
      size_t numUsed = 0;
      uint8_t idxInSib = 0;
      auto borderNode = root_->getLmBorderNode();
      while (reinterpret_cast<uintptr_t>(borderNode) != BTreeNodeT::NOTFOUND) {
        numUsed += reinterpret_cast<BtmNodeT *>(borderNode->getChildPtr(idxInSib))->getNumChildren();
        borderNode = getNextBtmRef(idxInSib);
      }
      return numUsed;
    }


    size_t calcNumSlotsBtmPart() const noexcept {
      size_t numUsed = 0;
      uint8_t idxInSib = 0;
      auto borderNode = root_->getLmBorderNode();
      while (reinterpret_cast<uintptr_t>(borderNode) != BTreeNodeT::NOTFOUND) {
        numUsed += kBtmB;
        borderNode = getNextBtmRef(idxInSib);
      }
      return numUsed;
    }


    void printStatictics(std::ostream & os) const noexcept {
      const size_t totalLen = getSumOfWeight();
      const size_t numRuns = calcNumRuns();
      os << "TotalLen = " << totalLen << ", #Runs = " << numRuns << ", Alphabet Size = " << calcNumAlph() << ", BTree arity param B = " << static_cast<int>(B) << std::endl;
      os << "Total: " << calcMemBytes() << " bytes" << std::endl;
      os << "MTree: " << calcMemBytesMTree() << " bytes, OccuRate = " << ((rootM_->calcNumSlots()) ? 100.0 * rootM_->calcNumUsed() / rootM_->calcNumSlots() : 0)
         << " (= 100*" << rootM_->calcNumUsed() << "/" << rootM_->calcNumSlots() << ")" << std::endl;
      os << "ATree: " << calcMemBytesATree() << " bytes, OccuRate = " << ((rootA_->calcNumSlots()) ? 100.0 * rootA_->calcNumUsed() / rootA_->calcNumSlots() : 0)
         << " (= 100*" << rootA_->calcNumUsed() << "/" << rootA_->calcNumSlots() << ")" << std::endl;
      os << "STree: " << calcMemBytesSTree() << " bytes, OccuRate = " << ((calcNumSlotsSTree()) ? 100.0 * calcNumUsedSTree() / calcNumSlotsSTree() : 0)
         << " (= 100*" << calcNumUsedSTree() << "/" << calcNumSlotsSTree() << ")" << std::endl;
      os << "IdxConvertVecs: " << calcMemBytesIdxConvertVecs() << " bytes ~ "
         << "(2*" << static_cast<int>(idxM2S_.getW()) << "(bitwidth)*" << idxM2S_.capacity() << "(capacity each))/8, "
         << "OccuRate = " << ((idxM2S_.capacity() + idxS2M_.capacity()) ? 100.0 * 2 * numRuns / (idxM2S_.capacity() + idxS2M_.capacity()) : 0)
         << " (= 100*2*" << numRuns << "/" << (idxM2S_.capacity() + idxS2M_.capacity()) << ")" << std::endl;
      os << "WeightVecs: " << calcMemBytesWeightVecs() << " bytes" << std::endl;
      os << "BtmArrays: " << calcMemBytesBtmArrays() << " bytes, "
         << "OccuRate = " << ((idxM2S_.capacity() + idxS2M_.capacity()) ? 100.0 * (idxM2S_.size() + idxS2M_.size()) / (idxM2S_.capacity() + idxS2M_.capacity()) : 0)
         << " (= 100*" << (idxM2S_.size() + idxS2M_.size())/B << "/" << (idxM2S_.capacity() + idxS2M_.capacity())/B << "), "
         << "OccuRate (btmM) = " << ((idxM2S_.capacity()) ? 100.0 * idxM2S_.size() / idxM2S_.capacity() : 0)
         << " (= 100*" << idxM2S_.size()/B << "/" << idxM2S_.capacity()/B << "), "
         << "OccuRate (btmS) = " << ((idxS2M_.capacity()) ? 100.0 * idxS2M_.size() / idxS2M_.capacity() : 0)
         << " (= 100*" << idxS2M_.size()/B << "/" << idxS2M_.capacity()/B << ")" << std::endl;
    }


    void printDebugInfo(std::ostream & os) const noexcept {
      {
        uint64_t c = UINT64_MAX;
        std::cout << "check runs:" << std::endl;
        uint64_t pos = 0;
        uint64_t len = 0;
        for (auto idxM = searchPosM(pos); idxM != BTreeNode<B>::NOTFOUND; idxM = getNextIdxM(idxM)) {
          ++pos;
          len += getWeightFromIdxM(idxM);
          if (getWeightFromIdxM(idxM) == 0) {
            std::cout << "detected 0 length run: " << idxM << ", " << pos << std::endl;
          }
          if (c == getCharFromIdxM(idxM)) {
            auto idxM0 = getPrevIdxM(idxM);
            std::cout << "detected consecutive runs having the same char: " 
                      << idxM << ", " << pos << ", (" << c << ", " << getWeightFromIdxM(idxM0) << ")" << ", (" << c << ", " << getWeightFromIdxM(idxM) << ")" << std::endl;
          }
          c = getCharFromIdxM(idxM);
        }
        std::cout << "run: " << pos << ", len: " << len << std::endl;
      }

      {
        uint64_t pos = 0;
        for (auto idxM = searchPosM(pos); idxM != BTreeNode<B>::NOTFOUND; idxM = getNextIdxM(idxM)) {
          os << "(" << idxM << ":" << getCharFromIdxM(idxM) << "^" << getWeightFromIdxM(idxM) << ") ";
        }
        os << std::endl;
      }

      {
        const uint64_t numBtmM = idxM2S_.size() / B;
        os << "information on M" << std::endl;
        for (uint64_t i = 0; i < numBtmM; ++i) {
          const auto nextBtmM = getNextBtmM(i);
          os << "[" << i*B << "-" << (i+1)*B-1 << "] (num=" << (int)getNumChildrenM(i) << " lbl=" 
             << labelM_[i] << " par=" << parentM_[i] << " sib=" << (int)idxInSiblingM_[i] << ") "
             << "=> " << nextBtmM * B << std::endl;
          for (uint64_t j = 0; j < getNumChildrenM(i); ++j) {
            if (j < getNumChildrenM(i) && B*i+j != idxS2M_.read(idxM2S_.read(B*i+j))) {
              os << "!!"; // WARNING, links are not maintained correctly
            }
            os << idxM2S_.read(B*i+j) << "(" << getWeightFromIdxM(B*i+j) << ")  ";
          }
          os << std::endl;
        }
      }

      {
        const uint64_t numBtmS = idxS2M_.size() / B;
        os << "information on S" << std::endl;
        for (uint64_t i = 0; i < numBtmS; ++i) {
          const auto nextIdxS = getNextIdxS(i*B + numChildrenS_[i] - 1);
          os << "[" << i*B << "-" << (i+1)*B-1 << "] (num=" << (int)numChildrenS_[i] << " ch=" << charS_[i] << " par=" 
             << parentS_[i] << " sib=" << (int)idxInSiblingS_[i] << ") "
             << "=> " << nextIdxS << std::endl;
          for (uint64_t j = 0; j < B; ++j) {
            os << idxS2M_.read(B*i+j) << "  ";
          }
          os << std::endl;
        }
      }

      os << "Alphabet: " << std::endl;
      for (const auto * rootS = getFstRootS();
           reinterpret_cast<uintptr_t>(rootS) != BTreeNode<B>::NOTFOUND;
           rootS = getNextRootS(rootS)) {
        const uint64_t btmS = reinterpret_cast<uintptr_t>(rootS->getLmBtm());
        os << "(" << charS_[btmS] << ", " << rootS->getSumOfWeight(kRow0) << ") ";
      }
      os << std::endl;
    }
  };
} // namespace itmmti

#endif
