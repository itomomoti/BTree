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
#ifndef INCLUDE_GUARD_PSumWithValue
#define INCLUDE_GUARD_PSumWithValue

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
  ////////////////////////////////////////////////////////////////
  class PackedArray_Null
  {
  public:
    using ValT = uint64_t;


  public:
    PackedArray_Null()
    {}


    ~PackedArray_Null()
    {}


    void clearAll() {
    }


    uint64_t getW() const noexcept {
      return 0;
    }


    void resize
    (
     size_t newSize
     ) noexcept {
    }


    /*!
     * @brief write sample
     */
    uint64_t read
    (
     const uint64_t idx
     ) const noexcept {
      return 0;
    }


    /*!
     * @brief write sample
     */
    void write
    (
     const uint64_t newSample,
     const uint64_t idx
     ) noexcept {
    }


    void increaseW
    (
     const uint64_t givenUpperBound
     ) {
    }


    void decreaseW
    (
     const uint64_t givenUpperBound
     ) {
    }


    void mvVals
    (
     const uint64_t srcIdx,
     const uint64_t tgtIdx,
     const uint8_t num
     ) noexcept {
    }


    //////////////////////////////// statistics
    size_t calcMemBytes
    (
     bool includeThis = true
     ) const noexcept {
      return 0;
    }


    size_t calcMemBytes_used() const noexcept {
      return 0;
    }
  };




  ////////////////////////////////////////////////////////////////
  template<uint16_t kBtmB>
  class BtmNodeForPSumWithVal
  {
  public:
    //// Public constant, alias etc.
    using BtmNodeT = BtmNodeForPSumWithVal<kBtmB>;


  private:
    //// Private member variables.
    StepCodeCore<2 * kBtmB> stcc_; //!< weights and values are stored alternatingly.
    uint16_t bitCapacity_; //!< Current bit capacity.
    uint16_t bitSize_; //!< Current bit size.
    uint8_t numChildren_; //!< Current size (number of elements).
    uint8_t wCodesAuxM_[2 * kBtmB / StepCodeUtil::kWCNum - 1];


  public:
    BtmNodeForPSumWithVal
    (
     uint16_t initBitCapacity = 0
     ) : stcc_(), bitCapacity_(0), bitSize_(0), numChildren_(0)
    {
      assert(initBitCapacity <= UINT16_MAX - 64);

      bitCapacity_ = static_cast<uint16_t>(stcc_.setBitCapacity(initBitCapacity));
    }


    ~BtmNodeForPSumWithVal()
    {
      // stcc_ is freed.
    }


    uint16_t getBitSize() const noexcept
    {
      return bitSize_;
    }


    uint16_t getBitCapacity() const noexcept
    {
      return bitCapacity_;
    }


    /*!
     * @brief Change capacity to max of givenCapacity and current size.
     * @node If givenCapacity is 0, it works as shrink_to_fit.
     */
    uint16_t changeBitCapacity
    (
     uint16_t givenBitCapacity
     ) {
      if (bitCapacity_ != givenBitCapacity) {
        bitCapacity_ = static_cast<uint16_t>(stcc_.setBitCapacity(std::max(bitSize_, givenBitCapacity)));
      }
    }


    uint16_t getNumChildren() const noexcept
    {
      return numChildren_;
    }


    uint64_t getWeight
    (
     uint16_t childIdx
     ) const noexcept {
      assert(childIdx < static_cast<uint16_t>(numChildren_));

      return stcc_.read(2 * childIdx);
    }


    uint64_t getVal
    (
     uint16_t childIdx
     ) const noexcept {
      assert(childIdx < static_cast<uint16_t>(numChildren_));

      return stcc_.read(2 * childIdx + 1);
    }


    uint64_t calcSumOfWeight
    (
     uint16_t childIdx_beg,
     uint16_t childIdx_end
     ) const noexcept {
      assert(childIdx_beg < childIdx_end);
      assert(childIdx_end <= static_cast<uint16_t>(numChildren_));

      uint64_t sum = 0;
      uint64_t bitPos = calcBitPos(2 * childIdx_beg);
      for (uint16_t i = childIdx_beg; i < childIdx_end; ++i) {
        const auto w = stcc_.readW(2 * i);
        sum += stcc_.readWBits(bitPos, w);
        bitPos += w + stcc_.readW(2 * i + 1);
      }

      return sum;
    }


    uint64_t calcSumOfWeight() const noexcept
    {
      uint64_t sum = 0;
      uint64_t bitPos = 0;
      for (uint64_t i = 0; i < numChildren_; ++i) {
        const auto w = stcc_.readW(2 * i);
        sum += stcc_.readWBits(bitPos, w);
        bitPos += w + stcc_.readW(2 * i + 1);
      }

      return sum;
    }


    /*!
     * @brief Calculate the beginning bit-pos of "idx"-th value in stcc_.
     */
    uint16_t calcBitPos
    (
     const uint16_t idx //!< in [0..2*numChildren_]
     ) const noexcept {
      assert(idx <= 2 * static_cast<uint16_t>(numChildren_));

      if (idx < 2 * static_cast<uint16_t>(numChildren_)) {
        return static_cast<uint16_t>(stcc_.calcBitPos(idx, wCodesAuxM_));
      } else {
        return bitSize_;
      }
    }


    /*!
     * @brief Return child index where partial sum of "pos+1" is achieved.
     * @pre Answer should be found in this bottom node.
     */
    uint8_t searchPos
    (
     uint64_t & pos, //!< [in,out] Give position to search. It is modified to relative position.
     uint64_t & retWeight, //!< [out] Capture last weight where pos exists.
     uint64_t & retVal //!< [out] Capture value.
     ) const noexcept {
      assert(pos < calcSumOfWeight()); // pre

      uint16_t bitPos = 0;
      uint16_t i = 0;
      while (true) {
        const auto w = stcc_.readW(i);
        retWeight = stcc_.readWBits(bitPos, w);
        const auto w_val = stcc_.readW(++i);
        if (pos < retWeight) {
          retVal = stcc_.readWBits(bitPos + w, w_val);
          return static_cast<uint8_t>(i / 2);
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
    uint8_t searchPos
    (
     uint64_t & pos, //!< [in,out] Give position to search. It is modified to relative position.
     uint64_t & retWeight //!< [out] Capture last weight where pos exists.
     ) const noexcept {
      assert(pos < calcSumOfWeight()); // pre

      uint16_t bitPos = 0;
      uint16_t i = 0;
      while (true) {
        const auto w = stcc_.readW(i);
        retWeight = stcc_.readWBits(bitPos, w);
        if (pos < retWeight) {
          return static_cast<uint8_t>(i / 2);
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
    uint8_t searchPos
    (
     uint64_t & pos //!< [in,out] Give position to search. It is modified to relative position.
     ) const noexcept {
      assert(pos < calcSumOfWeight()); // pre

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


    uint16_t setBitCapacity
    (
     uint16_t givenBitCapacity
     ) {
      return static_cast<uint16_t>(stcc_.setBitCapacity(static_cast<size_t>(givenBitCapacity)));
    }


    /*!
     * @brief Resize "numChildren_" to "newSize".
     * @note
     *   It does not change bitCapacity.
     */
    void resize
    (
     const uint8_t newSize
     ) noexcept {
      assert(newSize <= kBtmB);

      numChildren_ = newSize;
    }


    void changeWCodes
    (
     const uint64_t * srcWCodes,
     const uint16_t srcIdxBeg, //!< Beginning idx of src.
     const uint16_t srcLen, //!< Length of wCodes of src to insert.
     const uint16_t tgtIdxBeg, //!< Beginning idx of tgt.
     const uint16_t tgtLen //!< Length of wCodes of tgt to delete.
     ) noexcept {
      assert(tgtIdxBeg + tgtLen <= static_cast<uint16_t>(numChildren_) * 2);
      assert(tgtIdxBeg - tgtLen + srcLen <= kBtmB * 2);

      uint16_t sizeTimes2 = static_cast<uint16_t>(numChildren_) * 2;
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
     const uint16_t idxEnd
     ) noexcept {
      assert(idxBeg < idxEnd);

      const uint64_t beg = idxBeg / StepCodeUtil::kWCNum;
      const uint64_t end = (idxEnd - 1) / StepCodeUtil::kWCNum + (idxEnd <= (2 * kBtmB - StepCodeUtil::kWCNum));
      // std::cout << __FUNCTION__ << ": " << idxBeg << "->" << beg << ", " << idxEnd << "->" << end << std::endl;
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
      assert(childIdx + numChild_del <= numChildren_);
      {//debug
        std::cout << __FUNCTION__ << std::endl;
      }

      std::tuple<uint64_t, uint64_t, uint64_t> changeList[2];
      uint8_t clSize = 0;
      const uint16_t numL_old = lnode->getNumChildren();
      const uint16_t numR_old = this->getNumChildren();
      const uint16_t numTotal = numL_old + numR_old + numChild_ins - numChild_del;
      const uint16_t numL_new = numTotal / 2;
      const uint16_t numR_new = numTotal - numL_new;
      const uint16_t numToLeft = numL_new - numL_old;
      const bool isNewElemInL = childIdx < numToLeft;
      const bool isNewElemInR = childIdx + numChild_ins > numToLeft;
      uint16_t sumWL = lnode->bitSize_;
      uint16_t numL = numL_old;
      uint16_t curNumAfterDel = 0;
      uint16_t curNumSrcWCodes = 0;
      {
        const auto num = (isNewElemInL)? childIdx : numToLeft;
        if (num) {
          const auto w = this->calcBitPos(2 * num);
          changeList[clSize++] = {0, sumWL, w};
          sumWL += w;
          lnode->stcc_.mvWCodes(this->getConstPtr_wCodes(), 0, 2 * numL, 2 * num);
          numL += num;
        }
      }
      if (isNewElemInL) {
        curNumSrcWCodes = std::min(static_cast<uint16_t>(numToLeft - childIdx), numChild_ins);
        sumWL += StepCodeUtil::sumW(srcWCodes, 0, 2 * curNumSrcWCodes);
        lnode->stcc_.mvWCodes(srcWCodes, 0, 2 * numL, 2 * curNumSrcWCodes);
        numL += curNumSrcWCodes;
        if (numL < numL_new) { // Still need to move elements to left after inserting srcWCodes
          curNumAfterDel = numL_new - numL;
          const auto w = this->stcc_.sumW(2 * (childIdx + numChild_del), 2 * (childIdx + numChild_del + curNumAfterDel));
          changeList[clSize++] = {this->calcBitPos(2 * (childIdx + numChild_del)), sumWL, w};
          sumWL += w;
          lnode->stcc_.mvWCodes(this->getConstPtr_wCodes(), 2 * (childIdx + numChild_del), 2 * numL, 2 * curNumAfterDel);
        }
      }
      lnode->numChildren_ = static_cast<uint8_t>(numL_new);
      lnode->updateWCodesAuxM(2 * numL_old, 2 * numL_new);
      { // Update vals of lnode.
        if (lnode->bitCapacity_ < sumWL) {
          lnode->bitCapacity_ = lnode->setBitCapacity(sumWL);
        }
        for (uint8_t i = 0; i < clSize; ++i) {
          lnode->mvVals(this->getConstPtr_vals(), std::get<0>(changeList[i]), std::get<1>(changeList[i]), std::get<2>(changeList[i]));
        }
        lnode->bitSize_ = sumWL;
      }

      // Update this node.
      clSize = 0;
      const uint16_t bitPosOfLastChunk = this->calcBitPos(2 * (childIdx + numChild_del + curNumAfterDel));
      const uint16_t bitSizeOfLastChunk = this->bitSize_ - bitPosOfLastChunk;
      uint16_t sumWR = bitSizeOfLastChunk;
      if (numToLeft < childIdx) {
        const uint16_t num = childIdx - numToLeft;
        const uint16_t bitPos = this->calcBitPos(2 * numToLeft);
        const uint16_t w = this->calcBitPos(2 * childIdx) - bitPos;
        sumWR += w;
        changeList[clSize++] = {bitPos, 0, w};
        this->stcc_.mvWCodes(this->getConstPtr_wCodes(), 2 * numToLeft, 0, 2 * num);
      }
      if (isNewElemInR) {
        sumWR += StepCodeUtil::sumW(srcWCodes, 2 * curNumSrcWCodes, 2 * numChild_ins);
      }
      if (numR_old != childIdx + numChild_del) { // There are remaining children in tail.
        if (numR_old != numR_new) { // Need shift wCodes of "this" node.
          const uint16_t srcBeg = childIdx + numChild_del + curNumAfterDel;
          const uint16_t num = numR_old - srcBeg;
          const uint16_t tgtBeg = numR_new - num;
          this->stcc_.mvWCodes(this->getConstPtr_wCodes(), 2 * srcBeg, 2 * tgtBeg, 2 * num);
        }
        changeList[clSize++] = {bitPosOfLastChunk, sumWR - bitSizeOfLastChunk, bitSizeOfLastChunk};
      }
      if (isNewElemInR) {
        const uint16_t num = numChild_ins - curNumSrcWCodes;
        this->stcc_.mvWCodes(srcWCodes, 2 * curNumSrcWCodes, 2 * (childIdx + curNumSrcWCodes - numToLeft), 2 * num);
      }
      this->numChildren_ = static_cast<uint8_t>(numR_new);
      this->updateWCodesAuxM(0, 2 * numR_new);
      { // Update vals of "this" node.
        if (this->bitCapacity_ < sumWR) {
          this->bitCapacity_ = this->setBitCapacity(sumWR);
        }
        for (uint8_t i = 0; i < clSize; ++i) {
          std::cout << "mvVals: " << this->getConstPtr_vals() << ", " << std::get<0>(changeList[i]) << ", " << std::get<1>(changeList[i]) << ", " << std::get<2>(changeList[i]) << std::endl;
          this->mvVals(this->getConstPtr_vals(), std::get<0>(changeList[i]), std::get<1>(changeList[i]), std::get<2>(changeList[i]));
        }
        this->bitSize_ = sumWR;
      }
    }


    // void overflowToL_wCodes
    // (
    //  BtmNodeT * lnode,
    //  const uint64_t * srcWCodes,
    //  const uint16_t numChild_ins,
    //  const uint16_t childIdx, //!< Beginning childIdx of tgt.
    //  const uint16_t numChild_del //!< Length of wCodes of tgt to delete.
    //  ) noexcept {
    //   assert(childIdx + numChild_del <= numChildren_);

    //   std::tuple<uint64_t, uint64_t, uint64_t> changeList[2];
    //   uint8_t clSize = 0;
    //   const uint16_t numL_old = lnode->getNumChildren();
    //   const uint16_t numR_old = this->getNumChildren();
    //   const uint16_t numTotal = numL_old + numR_old + numChild_ins - numChild_del;
    //   const uint16_t numL_new = numTotal / 2;
    //   const uint16_t numR_new = numTotal - numL_new;
    //   const uint16_t numToLeft = numL_new - numL_old;

    //   uint16_t numSrcWCodesInL = 0;
    //   uint16_t numToLeft1 = 0; // Num of elements that exist before inserting position and to be moved to left node
    //   uint16_t numToLeft2 = 0; // Num of elements that exist after inserting position and to be moved to left node
    //   uint16_t numTailInR = numChildren_ - (childIdx + numChild_del); // Num of elements that exist after inserting position and to be remained in right node
    //   if (childIdx < numToLeft) { // new elements are in L
    //     numToLeft1 = childIdx;
    //     if (childIdx + numChild_ins <= numToLeft) { // new elements are only in L
    //       numSrcWCodesInL = numChild_ins;
    //       numToLeft2 = numToLeft - (childIdx + numChild_ins);
    //       numTailInR -= numToLeft2;
    //     } else { // new elements are also in in R
    //       numSrcWCodesInL = numToLeft - childIdx;
    //     }
    //   } else { // new elements are only in R
    //     numToLeft1 = childIdx - numL_new;
    //     numToLeft2 = numL_old - (childIdx + numChild_del);
    //   }

    //   uint16_t numL = lnode->getNumChildren();
    //   uint64_t sumWL = lnode->bitSize_;
    //   if (numToLeft1) {
    //     const auto bitPos = this->calcBitPos(2 * childIdx);
    //     const auto w = this->calcBitPos(2 * (childIdx + numToLeft1)) - bitPos;
    //     changeList[clSize++] = {bitPos, sumWL, w};
    //     sumWL += w;
    //     lnode->stcc_.mvWCodes(this->getConstPtr_wCodes(), 0, 2 * numL, 2 * numToLeft1);
    //     numL += numToLeft1;
    //   }
    //   if (numSrcWCodesInL) {
    //     sumWL += StepCodeUtil::sumW(srcWCodes, 0, 2 * numSrcWCodesInL);
    //     lnode->stcc_.mvWCodes(srcWCodes, 0, 2 * numSrcWCodesInL, 2 * numL, 2 * numSrcWCodesInL);
    //     numL += numSrcWCodesInL;
    //   }
    //   if (numToLeft2) {
    //     const auto bitPos = this->calcBitPos(2 * (childIdx + numChild_del));
    //     const auto w = this->bitSize_ - bitPos;
    //     changeList[clSize++] = {bitPos, sumWL, w};
    //     sumWL += w;
    //     lnode->stcc_.mvWCodes(this->getConstPtr_wCodes(), 2 * (childIdx + numChild_del), 2 * numL, 2 * numToLeft2);
    //   }
    //   lnode->numChildren_ = numL_new;
    //   lnode->updateWCodesAuxM(2 * numL_old, 2 * numL_new);
    //   { // Update vals of left node
    //     if (lnode->bitCapacity_ < sumWL) {
    //       lnode->bitCapacity_ = lnode->setBitCapacity(sumWL);
    //     }
    //     for (uint8_t i = 0; i < clSize; ++i) {
    //       lnode->mvVals(this->getConstPtr_vals(), std::get<0>(changeList[i]), std::get<1>(changeList[i]), std::get<2>(changeList[i]));
    //     }
    //     lnode->bitSize_ += sumWL;
    //   }

    //   if (numTailInR) {
        
    //   } else {
        
    //   }

    //   if (numSrcWCodesInL) {
    //     const uint64_t sumWL_ins = StepCodeUtil::sumW(srcWCodes, 0, 2 * numSrcWCodesInL);
    //     const uint64_t tailBitPos_new = this->calcBitPos(2 * childIdx) + sumWL_ins;
    //     this->bitSize_ = tailBitPos_new;
    //     const auto numTail = numL_new - (childIdx + numChild_ins);
    //     if (numTail) {
    //       const auto tailBitPos_old = this->calcBitPos(2 * (childIdx + numChild_del));
    //       const auto w = this->calcBitPos(2 * (childIdx + numChild_del + numTail)) - tailBitPos_old;
    //       this->bitSize_ += w;
    //       if (tailBitPos_new != tailBitPos_old) {
    //         if (this->bitCapacity_ < this->bitSize_) {
    //           this->bitCapacity_ = this->setBitCapacity(this->bitSize_);
    //         }
    //         this->mvVals(this->getConstPtr_vals(), tailBitPos_old, tailBitPos_new, w);
    //       }
    //       if (numChild_ins != numChild_del) {
    //         this->stcc_.mvWCodes(this->getConstPtr_wCodes(), 2 * (childIdx * numChild_del), 2 * (childIdx * numChild_ins), 2 * numTail);
    //       }
    //     } else {
    //       if (this->bitCapacity_ < this->bitSize_) {
    //         this->bitCapacity_ = this->setBitCapacity(this->bitSize_);
    //       }
    //     }
    //     this->stcc_.mvWCodes(srcWCodes, 0, 2 * childIdx, 2 * numSrcWCodesInL);
    //     this->updateWCodesAuxM(2 * childIdx, 2 * numL_new);
    //   } else {
    //     this->bitSize_ = this->calcBitPos(2 * numL_new); // shrink (just change bitSize)
    //   }
    //   this->numChildren_ = numL_new;
    // }


    void overflowToR_wCodes
    (
     BtmNodeT * rnode,
     const uint64_t * srcWCodes,
     const uint16_t numChild_ins,
     const uint16_t childIdx, //!< Beginning childIdx of tgt.
     const uint16_t numChild_del //!< Length of wCodes of tgt to delete.
     ) noexcept {
      assert(childIdx + numChild_del <= numChildren_);
      {//debug
        std::cout << __FUNCTION__ << std::endl;
      }

      std::tuple<uint64_t, uint64_t, uint64_t> changeList[2];
      uint8_t clSize = 0;
      const uint16_t numL_old = this->getNumChildren();
      const uint16_t numR_old = rnode->getNumChildren();
      const uint16_t numTotal = numL_old + numR_old + numChild_ins - numChild_del;
      const uint16_t numL_new = numTotal / 2;
      const uint16_t numR_new = numTotal - numL_new;
      const uint16_t numToRight = numR_new - numR_old;

      uint16_t numSrcWCodesInL = 0;
      uint16_t numToRight1 = 0;
      uint16_t numToRight2 = 0;
      if (childIdx < numL_new) { // new elements are in L
        if (childIdx + numChild_ins <= numL_new) { // new elements are only in L
          numSrcWCodesInL = numChild_ins;
          numToRight2 = numToRight;
        } else { // new elements are also in in R
          numSrcWCodesInL = numL_new - childIdx;
          numToRight2 = numL_old - (childIdx + numChild_del);
          // {
          //   std::cout << "koko: numToRight2 = " << numToRight2 << std::endl;
          // }
        }
      } else { // new elements are in R
        numToRight1 = childIdx - numL_new;
        numToRight2 = numL_old - (childIdx + numChild_del);
      }

      if (numR_old) { // shift wCodes of R to make space
        rnode->stcc_.mvWCodes(rnode->getConstPtr_wCodes(), 0, 2 * numToRight, 2 * numR_old);
      }

      uint16_t numR_increment = 0;
      uint16_t sumWR_increment = 0;
      if (numToRight1) {
        const uint16_t bitPos = this->calcBitPos(2 * (childIdx - numToRight1));
        const uint16_t w = this->calcBitPos(2 * childIdx) - bitPos;
        changeList[clSize++] = {bitPos, 0, w};
        sumWR_increment += w;
        rnode->stcc_.mvWCodes(this->getConstPtr_wCodes(), 2 * (childIdx - numToRight1), 0, 2 * numToRight1);
        numR_increment += numToRight1;
      }
      if (numSrcWCodesInL != numChild_ins) {
        sumWR_increment += StepCodeUtil::sumW(srcWCodes, 2 * numSrcWCodesInL, 2 * numChild_ins);
        rnode->stcc_.mvWCodes(srcWCodes, 2 * numSrcWCodesInL, 2 * numR_increment, 2 * (numChild_ins - numSrcWCodesInL));
        numR_increment += (numChild_ins - numSrcWCodesInL);
      }
      if (numToRight2) {
        const uint16_t bitPos = this->calcBitPos(2 * (numL_old - numToRight2));
        const uint16_t w = this->bitSize_ - bitPos;
        changeList[clSize++] = {bitPos, sumWR_increment, w};
        sumWR_increment += w;
        rnode->stcc_.mvWCodes(this->getConstPtr_wCodes(), 2 * (numL_old - numToRight2), 2 * numR_increment, 2 * numToRight2);
      }
      rnode->numChildren_ = static_cast<uint8_t>(numR_new);
      rnode->updateWCodesAuxM(0, 2 * numR_new);
      { // Update vals of "rnode".
        if (rnode->bitCapacity_ < rnode->bitSize_ + sumWR_increment) {
          rnode->bitCapacity_ = rnode->setBitCapacity(rnode->bitSize_ + sumWR_increment);
        }
        if (numR_old) {
          rnode->mvVals(rnode->getConstPtr_vals(), 0, sumWR_increment, rnode->bitSize_);
        }
        for (uint8_t i = 0; i < clSize; ++i) {
          rnode->mvVals(this->getConstPtr_vals(), std::get<0>(changeList[i]), std::get<1>(changeList[i]), std::get<2>(changeList[i]));
        }
        rnode->bitSize_ += sumWR_increment;
      }

      if (numSrcWCodesInL) {
        // {
        //   std::cout << "numSrcWCodesInL = " << numSrcWCodesInL << std::endl;
        // }
        const uint16_t sumWL_ins = static_cast<uint16_t>(StepCodeUtil::sumW(srcWCodes, 0, 2 * numSrcWCodesInL));
        const uint16_t tailBitPos_new = this->calcBitPos(2 * childIdx) + sumWL_ins;
        this->bitSize_ = tailBitPos_new;
        const uint16_t numTail = numL_new - (childIdx + numSrcWCodesInL);
        if (numTail) {
          const uint16_t tailBitPos_old = this->calcBitPos(2 * (childIdx + numChild_del));
          const uint16_t w = this->calcBitPos(2 * (childIdx + numChild_del + numTail)) - tailBitPos_old;
          this->bitSize_ += w;
          if (tailBitPos_new != tailBitPos_old) {
            if (this->bitCapacity_ < this->bitSize_) {
              this->bitCapacity_ = this->setBitCapacity(this->bitSize_);
            }
            this->mvVals(this->getConstPtr_vals(), tailBitPos_old, tailBitPos_new, w);
          }
          if (numChild_ins != numChild_del) {
            this->stcc_.mvWCodes(this->getConstPtr_wCodes(), 2 * (childIdx + numChild_del), 2 * (childIdx + numChild_ins), 2 * numTail);
          }
        } else {
          if (this->bitCapacity_ < this->bitSize_) {
            this->bitCapacity_ = this->setBitCapacity(this->bitSize_);
          }
        }
        this->stcc_.mvWCodes(srcWCodes, 0, 2 * childIdx, 2 * numSrcWCodesInL);
        this->updateWCodesAuxM(2 * childIdx, 2 * numL_new);
      } else {
        this->bitSize_ = this->calcBitPos(2 * numL_new); // shrink (just change bitSize)
        this->updateWCodesAuxM(2 * numL_new - 1, 2 * numL_new);
      }
      this->numChildren_ = static_cast<uint8_t>(numL_new);
    }


    void makeSpace
    (
     const uint16_t childIdx, //!< Beginning childIdx of tgt.
     const uint64_t * srcWCodes,
     const uint16_t sumW_ins,
     const uint16_t sumW_del,
     const uint8_t numChild_ins,
     const uint8_t numChild_del //!< Length of wCodes of tgt to delete.
     ) noexcept {
      // {//debug
      //   std::cerr << __FUNCTION__ << ": idxBase = " << idxBase << ", childIdx = " << (int)childIdx << ", sumW_ins = " << (int)sumW_ins
      //             << ", numChild_ins = " << (int)numChild_ins << ", numChild_del = " << (int)numChild_del << ", insertMorS = " << insertMorS << std::endl;
      //   // btmPtrs_[insertMorS][idxBase / kBtmB]->printDebugInfo(std::cerr);
      // }
      const uint16_t tailNum = this->numChildren_ - (childIdx + numChild_del); // at least 0 by assumption.
      uint16_t tailW = 0;
      if (tailNum) {
        tailW = this->stccSize_ - this->calcBitPos(2 * (childIdx + numChild_del));
        this->stcc_.mvWCodes(this->stcc_.getConstPtr_wCodes(), 2 * (childIdx + numChild_del), 2 * (childIdx + numChild_ins), 2 * tailNum);
      }
      this->stcc_.mvWCodes(srcWCodes, 0, 2 * childIdx, 2 * numChild_ins);
      this->numChildren_ += numChild_ins - numChild_del;
      this->updateWCodesAuxM(2 * childIdx, 2 * btmNodeM.numChildren_);
      if (sumW_ins != sumW_del) {
        const uint16_t newBitSize = this->stccSize_ + sumW_ins - sumW_del;
        this->reserveBitCapacity(newBitSize);
        if (tailNum) {
          this->stcc_.mvVals(this->stcc_.getConstPtr_vals(), this->stccSize_ - tailW, newBitSize - tailW, tailW);
        }
        this->stccSize_ = newBitSize;
      }
    }


    void writeNewElem
    (
     const uint16_t childIdx, //!< Beginning childIdx
     const uint64_t * newWeights, //!< Storing weights to insert to stcc
     const uint64_t * newVals, //!< Storing new values to insert
     const uint8_t numChild_ins
     ) noexcept {
      // {//debug
      //   std::cerr << __FUNCTION__ << ": lnode = " << lnode << ", rnode = " << rnode
      //             << ", childIdx = " << (int)childIdx << ", numChild_ins = " << (int)numChild_ins << std::endl;
      // }
      assert(childIdx + numChild_ins <= kBtmB);

      uint64_t bitPos = this->calcBitPos(2 * childIdx);
      for (uint16_t i = childIdx; i < childIdx + numChild_ins; ++i) {
        uint8_t w = this->stcc_.readW(2 * i);
        this->stcc_.writeWBits(newWeights[i - childIdx], bitPos, w);
        bitPos += w;
        w = this->stcc_.readW(2 * i + 1);
        this->stcc_.writeWBits(newVals[i - childIdx], bitPos, w);
        bitPos += w;
      }
    }


    void writeNewElemInTwo
    (
     BtmNodeT * rnode,
     const uint16_t childIdx_ins, //!< Beginning childIdx relative from beginning position of "this" node
     const uint64_t * newWeights, //!< Storing weights to insert to stcc
     const uint64_t * newVals, //!< Storing new values to insert
     const uint8_t numChild_ins
     ) {
      auto lnode = this;
      uint16_t numL = lnode->getNumChildren();
      uint16_t numNewElemInL = 0;
      if (childIdx_ins < numL) { // Insert new elements to lnode
        uint64_t bitPos = lnode->calcBitPos(2 * childIdx_ins);
        numNewElemInL = std::min(numChild_ins, static_cast<uint16_t>(numL - childIdx_ins));
        // {//debug
        //   std::cout << "insert to l: (" << lnode << ")" << numNewElemInL << std::endl;
        //   // lnode->printStatistics(std::cout, true);
        // }
        for (uint16_t ci = childIdx_ins; ci < childIdx_ins + numNewElemInL; ++ci) {
          uint8_t w = lnode->stcc_.readW(2 * ci);
          lnode->stcc_.writeWBits(newWeights[ci - childIdx_ins], bitPos, w);
          bitPos += w;
          w = lnode->stcc_.readW(2 * ci + 1);
          lnode->stcc_.writeWBits(newVals[ci - childIdx_ins], bitPos, w);
          bitPos += w;
        }
      }
      if (numNewElemInL < numChild_ins) { // Insert new elements to rnode
        // {//debug
        //   std::cout << "insert to r:" << std::endl;
        //   rnode->printStatistics(std::cout, true);
        // }
        childIdx_ins += numNewElemInL - numL;
        uint64_t bitPos = rnode->calcBitPos(2 * childIdx_ins);
        for (uint16_t ci = childIdx_ins; ci < childIdx_ins + numChild_ins - numNewElemInL; ++ci) {
          uint8_t w = rnode->stcc_.readW(2 * ci);
          // std::cout << "ci = " << ci
          //           << ", w = " << (int)w
          //           << ", weight = " << weights[ci - childIdx_ins] << std::endl;
          rnode->stcc_.writeWBits(newWeights[ci - childIdx_ins + numNewElemInL], bitPos, w);
          bitPos += w;
          w = rnode->stcc_.readW(2 * ci + 1);
          // std::cout << "ci = " << ci
          //           << ", w = " << (int)w
          //           << ", val = " << vals[ci - childIdx_ins] << std::endl;
          rnode->stcc_.writeWBits(newVals[ci - childIdx_ins + numNewElemInL], bitPos, w);
          bitPos += w;
        }
      }
    }


    /*!
     * @brief Merge elements in right node into this node
     * @note Right node is not deleted in this function
     */
    void merge
    (
     const BtmNodeT * rnode,
     const uint16_t childIdx, //!< Beginning childIdx relative from beginning position of "this" node
     const uint64_t * srcWCodes,
     const int16_t sumW_diff,
     const uint8_t numChild_ins,
     const uint8_t numChild_del //!< Length of wCodes of tgt to delete.
     ) noexcept {
      // {//debug
      //   std::cerr << __FUNCTION__ << ": idxBase = " << idxBase << ", childIdx = " << (int)childIdx << ", sumW_ins = " << (int)sumW_ins
      //             << ", numChild_ins = " << (int)numChild_ins << ", numChild_del = " << (int)numChild_del << ", insertMorS = " << insertMorS << std::endl;
      //   // btmPtrs_[insertMorS][idxBase / kBtmB]->printDebugInfo(std::cerr);
      // }
      uint16_t growStccSize = this->stccSize_;
      uint16_t growNum = this->numChildren_;
      {
        const uint16_t newBitSize = growStccSize + rnode->stccSize_ + sumW_diff;
        this->reserveBitCapacity(newBitSize);
        this->stccSize_ = newBitSize;
      }
      if (childIdx < growNum) {
        const uint16_t tailNum = growNum - (childIdx + numChild_del); // at least 0 by assumption.
        if (tailNum) {
          const uint16_t tailW = growStccSize - this->calcBitPos(2 * (childIdx + numChild_del));
          this->stcc_.mvWCodes(this->stcc_.getConstPtr_wCodes(), 2 * (childIdx + numChild_del), 2 * (childIdx + numChild_ins), 2 * tailNum);
          this->stcc_.mvVals(this->stcc_.getConstPtr_vals(), growStccSize - tailW, growStccSize + sumW_diff - tailW, tailW);
        }
        this->stcc_.mvWCodes(srcWCodes, 0, 2 * childIdx, 2 * numChild_ins);
        growNum += numChild_ins - numChild_del;
        growStccSize += sumW_diff;
        this->stcc_.mvWCodes(rnode->stcc_.getConstPtr_wCodes(), 0, 2 * growNum, 2 * rnode->numChildren_);
        this->stcc_.mvVals(rnode->stcc_.getConstPtr_vals(), 0, growStccSize, rnode->stccSize_);
        this->numChildren_ = growNum + rnode->numChildren_;
      } else {
        const uint16_t childIdxInR = childIdx - growNum;
        if (childIdxInR) {
          const uint16_t stccSize1 = rnode->calcBitPos(2 * childIdxInR1);
          this->stcc_.mvWCodes(rnode->stcc_.getConstPtr_wCodes(), 0, 2 * growNum, 2 * childIdxInR);
          this->stcc_.mvVals(rnode->stcc_.getConstPtr_vals(), 0, growStccSize, stccSize1);
          growNum += childIdxInR;
          growStccSize += stccSize1;
        }
        this->stcc_.mvWCodes(srcWCodes, 0, 2 * growNum, 2 * numChild_ins);
        growNum += numChild_ins;
        growStccSize += sumW_diff;
        const uint16_t tailNum = rnode->numChildren_ - (childIdxInR + numChild_del);
        if (tailNum) {
          const uint16_t tailW = rnode->stccSize_ - rnode->calcBitPos(2 * (childIdxInR + numChild_del));
          this->stcc_.mvWCodes(rnode->stcc_.getConstPtr_wCodes(), 2 * (childIdxInR + numChild_del), 2 * growNum, 2 * tailNum);
          this->stcc_.mvVals(rnode->stcc_.getConstPtr_vals(), rnode->stccSize_ - tailW, growStccSize, tailW);
        }
        this->numChildren_ = growNum + tailNum;
      }
    }


  public:
    template<class BTreeNodeT>
    void insert
    (
     const uint64_t * weights,
     const uint64_t * vals,
     const uint16_t numChild_ins,
     const uint16_t childIdx, //!< Beginning childIdx of tgt.
     const uint16_t numChild_del, //!< Length of wCodes of tgt to delete.
     BTreeNodeT * parent,
     const uint8_t idxInSibling
     ) {
      // {//debug
      //   std::cout << __func__ << ": numChild_ins = " << numChild_ins << ", childIdx = " << childIdx << ", numChild_del = " << numChild_del << std::endl;
      // }
      assert(numChild_ins <= kBtmB);
      assert(childIdx + numChild_del <= numChildren_);

      uint64_t wCodesTemp[2 * kBtmB / StepCodeUtil::kWCNum];
      uint16_t sumW_ins = 0;
      uint64_t sumWeights_ins = 0;
      for (uint16_t ci = 0; ci < numChild_ins; ++ci) {
        sumWeights_ins += weights[ci];
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
        uint64_t bitPos = this->calcBitPos(2 * childIdx);
        const uint64_t begPos = bitPos;
        for (uint16_t ci = childIdx; ci < childIdx + numChild_del; ++ci) {
          const uint8_t w = stcc_.readW(2 * ci);
          sumWeights_del += stcc_.readWBits(bitPos, w);
          bitPos += w + stcc_.readW(2 * ci + 1);
        }
        sumW_del = static_cast<uint16_t>(bitPos - begPos);
      }

      if (sumWeights_ins != sumWeights_del) {
        parent->changePSumFrom(idxInSibling, static_cast<int64_t>(sumWeights_ins) - static_cast<int64_t>(sumWeights_del));
      }

      if (numChildren_ - numChild_del + numChild_ins <= kBtmB) { // Easy case: This node can accommodate inserting elements.
        {//debug
          std::cout << __FUNCTION__ << ": easy case" << std::endl;
        }
        makeSpace(childIdx, srcWCodes, sumW_ins, sumW_del, numChild_ins, numChild_del);
        this->writeNewElem(childIdx, weights, vals, numChild_ins);
        return;
      }

      if (idxInSibling) { // Check previous sibling.
        lnode = reinterpret_cast<BtmNodeT *>(parent->getChildPtr(idxInSibling - 1));
        uint16_t num_old = lnode->getNumChildren();
        if (num_old <= kBtmB + numChild_del - numChild_ins) { // Previous sibling can accommodate overflowed elements.
          this->overflowToL_wCodes(lnode, wCodesTemp, numChild_ins, childIdx, numChild_del);
          lnode->writeNewElemInTwo(this, childIdx + num_old, weights, vals, numChild_ins);
          parent->changePSumAt(idxInSibling - 1, parent->getPSum(idxInSibling) + lnode->calcSumOfWeight(num_old, lnode->getNumChildren()));
          return;
        }
      }
      if (idxInSibling + 1 < parent->getNumChildren()) { // Check next sibling.
        auto rnode = reinterpret_cast<BtmNodeT *>(parent->getChildPtr(idxInSibling + 1));
        uint16_t num_old = rnode->getNumChildren();
        if (num_old <= kBtmB + numChild_del - numChild_ins) { // Next sibling can accommodate overflowed elements.
          this->overflowToR_wCodes(rnode, wCodesTemp, numChild_ins, childIdx, numChild_del);
          this->writeNewElemInTwo(rnode, childIdx, weights, vals, numChild_ins);
          parent->changePSumAt(idxInSibling, parent->getPSum(idxInSibling + 1) - rnode->calcSumOfWeight(0, rnode->getNumChildren() - num_old));
          return;
        }
      }

      { // This node has to be split
        {//debug
          std::cout << __FUNCTION__ << ": split" << std::endl;
        }
        auto rnode = new BtmNodeT();
        this->overflowToR_wCodes(rnode, wCodesTemp, numChild_ins, childIdx, numChild_del);
        this->writeNewElemInTwo(rnode, childIdx, weights, vals, numChild_ins);
        parent->handleSplitOfBtm(reinterpret_cast<BTreeNodeT *>(rnode), rnode->calcSumOfWeight(), idxInSibling);
        return;
      }
    }


    template<class BTreeNodeT>
    void insert_shrink
    (
     const uint64_t * weights,
     const uint64_t * vals,
     const uint16_t numChild_ins,
     const uint16_t childIdx, //!< Beginning childIdx of tgt.
     const uint16_t numChild_del, //!< Length of wCodes of tgt to delete.
     BTreeNodeT * parent,
     const uint8_t idxInSibling
     ) noexcept {
      // {//debug
      //   std::cout << __func__ << ": numChild_ins = " << numChild_ins << ", childIdx = " << childIdx << ", numChild_del = " << numChild_del << std::endl;
      // }
      assert(numChild_ins <= numChild_del);
      assert(childIdx + numChild_del <= numChildren_);

      uint64_t wCodesTemp[2 * kBtmB / StepCodeUtil::kWCNum];
      uint16_t sumW_ins = 0;
      uint64_t sumWeights_ins = 0;
      for (uint16_t ci = 0; ci < numChild_ins; ++ci) {
        sumWeights_ins += weights[ci];
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
        uint64_t bitPos = this->calcBitPos(2 * childIdx);
        const uint64_t begPos = bitPos;
        for (uint16_t ci = childIdx; ci < childIdx + numChild_del; ++ci) {
          const uint8_t w = stcc_.readW(2 * ci);
          sumWeights_del += stcc_.readWBits(bitPos, w);
          bitPos += w + stcc_.readW(2 * ci + 1);
        }
        sumW_del = static_cast<uint16_t>(bitPos - begPos);
      }

      if (sumWeights_ins != sumWeights_del) {
        parent->changePSumFrom(idxInSibling, static_cast<int64_t>(sumWeights_ins) - static_cast<int64_t>(sumWeights_del));
      }

      const uint16_t num_new = numChildren_ - numChild_del + numChild_ins;
      if (num_new == 0) {
        delete this;
        parent->handleDeleteOfBtm(0);
        return;
      }
      if (num_new < kBtmB / 2) {
        if (idxInSibling) { // Check previous sibling if it can be merged
          auto lnode = reinterpret_cast<BtmNodeT *>(parent->getChildPtr(idxInSibling - 1));
          const uint16_t lnum = lnode->getNumChildren();
          if (lnum + num_new <= 3 * kBtmB / 4) {
            lnode->merge(this, lnum + childIdx, srcWCodes, sumW_ins - sumW_del, numChild_ins, numChild_del);
            if (numChild_ins) {
              lnode->writeNewElem(lnum + childIdx, weights, vals, numChild_ins);
            }
            parent->changePSumAt(idxInSibling - 1, parent->getPSum(idxInSibling) + parent->getWeightOfChild(idxInSibling));
            parent->changePSumAt(idxInSibling, 0);
            delete this;
            parent->handleDeleteOfChild(0);
            return;
          }
        }
        if (idxInSibling + 1 < parent->getNumChildren()) { // Check next sibling if it can be merged
          auto rnode = reinterpret_cast<BtmNodeT *>(parent->getChildPtr(idxInSibling + 1));
          const uint16_t rnum = rnode->getNumChildren();
          if (rnum + num_new <= 3 * kBtmB / 4) {
            this->merge(rnode, childIdx, srcWCodes, sumW_ins - sumW_del, numChild_ins, numChild_del);
            if (numChild_ins) {
              this->writeNewElem(childIdx, weights, vals, numChild_ins);
            }
            parent->changePSumAt(idxInSibling, parent->getPSum(idxInSibling + 1) + parent->getWeightOfChild(idxInSibling + 1));
            parent->changePSumAt(idxInSibling + 1, 0);
            delete rnode;
            // parent->handleDeleteBtmNode();
            return;
          }
        }
      }

      {//debug
        std::cout << __FUNCTION__ << ": easy case" << std::endl;
      }
      makeSpace(childIdx, srcWCodes, sumW_ins, sumW_del, numChild_ins, numChild_del);
      if (numChild_ins) {
        this->writeNewElem(childIdx, weights, vals, numChild_ins);
      }
    }


    /*!
     * @brief Replace weights and/or values.
     */
    template<class BTreeNodeT>
    void replace
    (
     const uint64_t * array, //!< Storing weights and vals alternatingly (starting with weight iff "idx = 0 mod 2").
     const uint16_t num, //!< Number of elements to replace.
     const uint16_t idx, //!< in [0..2*numChildren_]. Beginning idx of tgt.
     BTreeNodeT * parent,
     const uint8_t idxInSibling
     ) {
      assert(idx + num <= 2 * static_cast<uint16_t>(numChildren_));

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
        if (sumW_ins > sumW_del && bitCapacity_ < bitSize_ + sumW_ins - sumW_del) {
          this->bitCapacity_ = this->setBitCapacity(bitSize_ + sumW_ins - sumW_del);
        }
        this->changeValPos(bitPos0, sumW_ins, sumW_del);
      }
      bitPos = bitPos0;
      for (uint16_t i = idx; i < idx + num; ++i) {
        uint8_t w = stcc_.readW(i);
        stcc_.writeWBits(array[i - idx], bitPos, w);
        bitPos += w;
      }

      if (weights_diff != 0) {
        parent->changePSumFrom(idxInSibling, weights_diff);
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


    void printStatistics
    (
     std::ostream & os,
     const bool verbose = false
     ) const noexcept {
      os << "BtmNodeForPSumWithVal object (" << this << ") " << __func__ << "(" << verbose << ") BEGIN" << std::endl;
      os << "BTree arity for bottom node = " << static_cast<int>(kBtmB) << std::endl;
      const size_t sumWeights = calcSumOfWeight();
      os << "num of children = " << static_cast<uint64_t>(numChildren_)
         << ", sum of weights = " << sumWeights << std::endl;
      os << "bit size = " << bitSize_ << ", bit capacity = " << bitCapacity_ << std::endl;
      os << "Total: " << calcMemBytes() << " bytes" << std::endl;
      os << "Memory usage for dynamic arrays of step code: " << calcMemBytesDynArray() << " bytes" << std::endl;
      if (verbose) {
        {
          const uint16_t num = static_cast<uint16_t>(numChildren_) * 2;
          std::cout << "dump bit witdth stored in wCodes (" << stcc_.getConstPtr_wCodes() << ")" << std::endl;
          for (uint64_t i = 0; i < num; ++i) {
            std::cout << (uint64_t)(stcc_.readW(i)) << " ";
          }
          std::cout << std::endl;
        }
        {
          const uint16_t num = static_cast<uint16_t>(numChildren_) * 2;
          std::cout << "dump values" << std::endl;
          for (uint64_t i = 0; i < num; ++i) {
            std::cout << stcc_.read(i) << " ";
          }
          std::cout << std::endl;
        }
        {
          std::cout << "dump bits in vals_ (" << stcc_.getConstPtr_vals() << ")" << std::endl;
          for (uint64_t i = 0; i < (bitSize_ + 63) / 64; ++i) {
            std::cout << "(" << i << ")";
            for (uint64_t j = 0; j < 64; ++j) {
              std::cout << bits::readWBits_S(stcc_.getConstPtr_vals(), 64 * i + 63 - j, ctcbits::UINTW_MAX(1));
            }
            std::cout << " ";
          }
          std::cout << std::endl;
        }
      }
      os << "BtmNodeForPSumWithVal object (" << this << ") " << __func__ << "(" << verbose << ") END" << std::endl;
    }
  };




  /*!
   * @brief Dynamic partial sum data structure implemented by B+tree, where each leaf is associated with value.
   * @tparam kB Arity for internal node of B+tree, which should be in {32, 64, 128}.
   * @tparam kBtmB Arity for bottom node of B+tree, which should be in {32, 64, 128}.
   * @par Notation
   */
  template<typename tparam_BTreeNodeT, typename LeafValT, typename tparam_BtmNodeT>
  class PSumWithValue
  {
  public:
    // Public constant, alias etc.
    using BTreeNodeT = tparam_BTreeNodeT;
    using BtmNodeT = tparam_BtmNodeT;
    static constexpr uint8_t kB{BTreeNodeT::kB};
    static constexpr uint8_t kBtmB{BtmNodeT::kBtmB};


  private:
    // Private member variables.
    typename BTreeNodeT::SuperRootT sr_; //!< Super root of B+tree.


  public:
    PSumWithValue()
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
      sr_.setRoot(new BTreeNodeT(nullptr, true, true, false, true, false));
      sr_.root_->putFirstBtm(reinterpret_cast<BTreeNodeT *>(fstBtmNode), 0);
      // Insert sentinel.
      const uint64_t dummyArray[] = {0};
      fstBtmNode->insert(dummyArray, dummyArray, 1, 0, 0, sr_.root_, 0);
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
      auto borderNode = sr_.root_->getLmBorderNode_DirectJump();
      while (reinterpret_cast<uintptr_t>(borderNode) != BTreeNodeT::NOTFOUND) {
        delete reinterpret_cast<BtmNodeT *>(borderNode->getChildPtr(idxInSib));
        borderNode = borderNode->getNextBtmRef_DirectJump(idxInSib);
      }
      // Delete uppdar nodes.
      sr_.root_->clearUpperNodes();
      sr_.root_ = nullptr;
    }


    /*!
     * @brief Return if data structure is ready.
     */
    bool isReady() const noexcept {
      return (sr_.root_ != NULL);
    }


    /*!
     * @brief Return |T|.
     */
    size_t getSumOfWeight() const noexcept {
      assert(isReady());

      return sr_.root_->getSumOfWeight();
    }


    BTreeNodeT * getRoot() noexcept {
      return sr_.root_;
    }


    const BTreeNodeT * getConstRoot() noexcept {
      return sr_.root_;
    }


  public:
    //// Search functions
    /*!
     * @brief Search for bottom node achieving psum of "pos+1".
     * @attention "pos" is modified to be the relative position (0base) from beginning of bottom.
     */
    BtmNodeT * searchBtm
    (
     uint64_t & pos, //!< [in,out] Give position to search (< |T|). It is modified to relative position.
     BTreeNodeT *& retNode, //!< [out] To capture parent of returned bottom node.
     uint8_t & retIdx //!< [out] To capture sibling idx of returned bottom node.
     ) const noexcept {
      assert(isReady());
      assert(pos < sr_.root_->getSumOfWeight());

      sr_.root_->searchPos(pos, retNode, retIdx);

      return reinterpret_cast<BtmNodeT *>(retNode->getChildPtr(retIdx));
    }


    BtmNodeT * getLmBtm
    (
     BTreeNodeT *& retNode, //!< [out] To capture parent of returned bottom node.
     uint8_t & retIdx //!< [out] To capture sibling idx of returned bottom node.
     ) const noexcept {
      assert(isReady());

      retNode = sr_.root_->getLmBorderNode_DirectJump();
      retIdx = 0;

      return reinterpret_cast<BtmNodeT *>(retNode->getChildPtr(0));
    }


    BtmNodeT * getRmBtm
    (
     BTreeNodeT *& retNode, //!< [out] To capture parent of returned bottom node.
     uint8_t & retIdx //!< [out] To capture sibling idx of returned bottom node.
     ) const noexcept {
      assert(isReady());

      retNode = sr_.root_->getRmBorderNode();
      retIdx = retNode->getNumChildren() - 1;

      return reinterpret_cast<BtmNodeT *>(retNode->getChildPtr(retIdx));
    }


  public:
    //// statistics
    size_t calcMemBytesUpperPart() const noexcept {
      return sr_.root_->calcMemBytes();
    }


    size_t calcMemBytesBtmPart() const noexcept {
      size_t size = 0;
      uint8_t idxInSib = 0;
      auto borderNode = sr_.root_->getLmBorderNode_DirectJump();
      while (reinterpret_cast<uintptr_t>(borderNode) != BTreeNodeT::NOTFOUND) {
        size += reinterpret_cast<BtmNodeT *>(borderNode->getChildPtr(idxInSib))->calcMemBytes();
        borderNode = borderNode->getNextBtmRef_DirectJump(idxInSib);
      }
      return size;
    }


    size_t calcMemBytesDynArrayOfStepCode() const noexcept {
      size_t size = 0;
      uint8_t idxInSib = 0;
      auto borderNode = sr_.root_->getLmBorderNode_DirectJump();
      while (reinterpret_cast<uintptr_t>(borderNode) != BTreeNodeT::NOTFOUND) {
        size += reinterpret_cast<BtmNodeT *>(borderNode->getChildPtr(idxInSib))->calcMemBytesDynArray();
        borderNode = borderNode->getNextBtmRef_DirectJump(idxInSib);
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
      return sr_.root_->calcNumUsed();
    }


    size_t calcNumSlotsUpperPart() const noexcept {
      return sr_.root_->calcNumSlots();
    }


    size_t calcNumUsedBtmPart() const noexcept {
      size_t numUsed = 0;
      uint8_t idxInSib = 0;
      auto borderNode = sr_.root_->getLmBorderNode_DirectJump();
      while (reinterpret_cast<uintptr_t>(borderNode) != BTreeNodeT::NOTFOUND) {
        numUsed += reinterpret_cast<BtmNodeT *>(borderNode->getChildPtr(idxInSib))->getNumChildren();
        borderNode = borderNode->getNextBtmRef_DirectJump(idxInSib);
      }
      return numUsed;
    }


    size_t calcNumSlotsBtmPart() const noexcept {
      size_t numSlots = 0;
      uint8_t idxInSib = 0;
      auto borderNode = sr_.root_->getLmBorderNode_DirectJump();
      while (reinterpret_cast<uintptr_t>(borderNode) != BTreeNodeT::NOTFOUND) {
        numSlots += kBtmB;
        borderNode = borderNode->getNextBtmRef_DirectJump(idxInSib);
      }
      return numSlots;
    }


    void printStatistics
    (
     std::ostream & os,
     const bool verbose = false
     ) const noexcept {
      os << "PSumWithValue object (" << this << ") " << __func__ << "(" << verbose << ") BEGIN" << std::endl;
      if (isReady()) {
        const size_t sumWeights = getSumOfWeight();
        os << "BTree arity = " << static_cast<int>(kB)
           << ", BTree arity for bottom node = " << static_cast<int>(kBtmB) << std::endl;
        os << "Sum of weights = " << sumWeights << std::endl;
        os << "Total: " << calcMemBytes() << " bytes" << std::endl;
        os << "Upper part of B+tree: " << calcMemBytesUpperPart()
           << " bytes, OccuRate = " << ((sr_.root_->calcNumSlots()) ? 100.0 * sr_.root_->calcNumUsed() / sr_.root_->calcNumSlots() : 0)
           << " (= 100*" << sr_.root_->calcNumUsed() << "/" << sr_.root_->calcNumSlots() << ")" << std::endl;
        os << "Bottom part of B+tree: " << calcMemBytesBtmPart()
           << " bytes, OccuRate = " << ((calcNumSlotsBtmPart()) ? 100.0 * calcNumUsedBtmPart() / calcNumSlotsBtmPart() : 0)
           << " (= 100*" << calcNumUsedBtmPart() << "/" << calcNumSlotsBtmPart() << ")" << std::endl;
        os << "Memory usage for dynamic arrays of step code: " << calcMemBytesDynArrayOfStepCode() << " bytes" << std::endl;
      } else {
        os << "Data structure is empty (not ready)." << std::endl;
      }
      os << "PSumWithValue object (" << this << ") " << __func__ << "(" << verbose << ") END" << std::endl;
    }


    // void printDebugInfo(std::ostream & os) const noexcept {
    //   {
    //     uint64_t c = UINT64_MAX;
    //     std::cout << "check runs:" << std::endl;
    //     uint64_t pos = 0;
    //     uint64_t len = 0;
    //     for (auto idxM = searchPosM(pos); idxM != BTreeNode<B>::NOTFOUND; idxM = getNextIdxM(idxM)) {
    //       ++pos;
    //       len += getWeightFromIdxM(idxM);
    //       if (getWeightFromIdxM(idxM) == 0) {
    //         std::cout << "detected 0 length run: " << idxM << ", " << pos << std::endl;
    //       }
    //       if (c == getCharFromIdxM(idxM)) {
    //         auto idxM0 = getPrevIdxM(idxM);
    //         std::cout << "detected consecutive runs having the same char: " 
    //                   << idxM << ", " << pos << ", (" << c << ", " << getWeightFromIdxM(idxM0) << ")" << ", (" << c << ", " << getWeightFromIdxM(idxM) << ")" << std::endl;
    //       }
    //       c = getCharFromIdxM(idxM);
    //     }
    //     std::cout << "run: " << pos << ", len: " << len << std::endl;
    //   }

    //   {
    //     uint64_t pos = 0;
    //     for (auto idxM = searchPosM(pos); idxM != BTreeNode<B>::NOTFOUND; idxM = getNextIdxM(idxM)) {
    //       os << "(" << idxM << ":" << getCharFromIdxM(idxM) << "^" << getWeightFromIdxM(idxM) << ") ";
    //     }
    //     os << std::endl;
    //   }

    //   {
    //     const uint64_t numBtmM = idxM2S_.size() / B;
    //     os << "information on M" << std::endl;
    //     for (uint64_t i = 0; i < numBtmM; ++i) {
    //       const auto nextBtmM = getNextBtmM(i);
    //       os << "[" << i*B << "-" << (i+1)*B-1 << "] (num=" << (int)getNumChildrenM(i) << " lbl=" 
    //          << labelM_[i] << " par=" << parentM_[i] << " sib=" << (int)idxInSiblingM_[i] << ") "
    //          << "=> " << nextBtmM * B << std::endl;
    //       for (uint64_t j = 0; j < getNumChildrenM(i); ++j) {
    //         if (j < getNumChildrenM(i) && B*i+j != idxS2M_.read(idxM2S_.read(B*i+j))) {
    //           os << "!!"; // WARNING, links are not maintained correctly
    //         }
    //         os << idxM2S_.read(B*i+j) << "(" << getWeightFromIdxM(B*i+j) << ")  ";
    //       }
    //       os << std::endl;
    //     }
    //   }

    //   {
    //     const uint64_t numBtmS = idxS2M_.size() / B;
    //     os << "information on S" << std::endl;
    //     for (uint64_t i = 0; i < numBtmS; ++i) {
    //       const auto nextIdxS = getNextIdxS(i*B + numChildrenS_[i] - 1);
    //       os << "[" << i*B << "-" << (i+1)*B-1 << "] (num=" << (int)numChildrenS_[i] << " ch=" << charS_[i] << " par=" 
    //          << parentS_[i] << " sib=" << (int)idxInSiblingS_[i] << ") "
    //          << "=> " << nextIdxS << std::endl;
    //       for (uint64_t j = 0; j < B; ++j) {
    //         os << idxS2M_.read(B*i+j) << "  ";
    //       }
    //       os << std::endl;
    //     }
    //   }

    //   os << "Alphabet: " << std::endl;
    //   for (const auto * rootS = getFstRootS();
    //        reinterpret_cast<uintptr_t>(rootS) != BTreeNode<B>::NOTFOUND;
    //        rootS = getNextRootS(rootS)) {
    //     const uint64_t btmS = reinterpret_cast<uintptr_t>(rootS->getLmBtm());
    //     os << "(" << charS_[btmS] << ", " << rootS->getSumOfWeight(kRow0) << ") ";
    //   }
    //   os << std::endl;
    // }
  };
} // namespace itmmti

#endif
