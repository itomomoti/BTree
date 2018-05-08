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
  /*!
   * @brief Bottom node for dynamic partial sum data structure implemented by B+tree, where each leaf is associated with value.
   * @tparam tparam_kBtmB Arity for bottom node of B+tree, which should be in {16, 32, 64, 128}.
   */
  template<uint16_t tparam_kBtmB>
  class BtmNodeForPSumWithVal
  {
  public:
    //// Public constant, alias etc.
    static constexpr uint16_t kBtmB{tparam_kBtmB};
    static constexpr size_t kUnitBits{kBtmB * 4};
    using BtmNodeT = BtmNodeForPSumWithVal<kBtmB>;


  private:
    //// Private member variables.
    StepCodeCore<2 * kBtmB> stcc_; //!< weights and values are stored alternatingly.
    uint16_t stccCapacity_; //!< Current bit capacity of stcc_
    uint16_t stccSize_; //!< Current bit size of stcc_
    uint8_t numChildren_; //!< Current size (number of elements).
    uint8_t wCodesAuxM_[2 * kBtmB / StepCodeUtil::kWCNum - 1];


  public:
    BtmNodeForPSumWithVal
    (
     uint16_t initStccCapacity = 0
     ) : stcc_(), stccCapacity_(0), stccSize_(0), numChildren_(0)
    {
      assert(initStccCapacity <= UINT16_MAX - 64);

      reserveBitCapacity(initStccCapacity);
    }


    ~BtmNodeForPSumWithVal()
    {
      // stcc_ is freed.
    }


    uint16_t getBitSize() const noexcept
    {
      return stccSize_;
    }


    uint16_t getBitCapacity() const noexcept
    {
      return stccCapacity_;
    }


    void reserveBitCapacity
    (
     uint16_t givenBitCapacity
     ) {
      if (givenBitCapacity > this->stccCapacity_) {
        size_t newSize = (static_cast<size_t>(givenBitCapacity) / kUnitBits + 2) * kUnitBits;
        this->stccCapacity_ = static_cast<uint16_t>(this->stcc_.setBitCapacity(static_cast<size_t>(givenBitCapacity)));
      }
    }


    void shrinkBitCapacity() {
      if (this->stccCapacity_ - this->stccSize_ > kUnitBits) {
        this->stccCapacity_ = static_cast<uint16_t>(this->stcc_.setBitCapacity(static_cast<size_t>(this->stccSize_)));
      }
    }


    uint16_t getNumChildren() const noexcept
    {
      return numChildren_;
    }


    uint64_t readNext
    (
     uint16_t idx, //!< Element idx to read in stcc_.
     uint64_t & bitPos //!< [in,out] BitPos of idx element in stcc_. Modified to BitPos of next element.
     ) const noexcept {
      const auto w = stcc_.readW(idx);
      const auto val = stcc_.readWBits(bitPos, w);
      bitPos += w;

      return val;
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
        return stccSize_;
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
     uint64_t & retVal, //!< [out] Capture value.
     uint64_t & bitPos //!< [out] Capture bitPos, pointing to end of read element.
     ) const noexcept {
      assert(pos < calcSumOfWeight()); // pre

      bitPos = 0;
      uint16_t i = 0;
      while (true) {
        const auto w = stcc_.readW(i);
        retWeight = stcc_.readWBits(bitPos, w);
        bitPos += w;
        const auto w_val = stcc_.readW(++i);
        if (pos < retWeight) {
          retVal = stcc_.readWBits(bitPos, w_val);
          bitPos += w_val;
          return static_cast<uint8_t>(i / 2);
        }
        bitPos += w_val;
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


    /*!
     * @brief Resize "numChildren_" to "newSize".
     * @note
     *   It does not change stccCapacity.
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
      stcc_.changeValPos(stccSize_, bitPos, insBitLen, delBitLen);
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
      assert(idxBeg <= idxEnd);
      assert(0 < idxEnd);

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
      {//debug
        std::cerr << __func__ << ": numChild_ins = " << numChild_ins << ", childIdx = " << childIdx << ", numChild_del = " << numChild_del << std::endl;
      }
      assert(childIdx + numChild_del <= numChildren_);

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
      uint16_t sumWL = lnode->stccSize_;
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
        lnode->reserveBitCapacity(sumWL);
        for (uint8_t i = 0; i < clSize; ++i) {
          lnode->mvVals(this->getConstPtr_vals(), std::get<0>(changeList[i]), std::get<1>(changeList[i]), std::get<2>(changeList[i]));
        }
        lnode->stccSize_ = sumWL;
      }

      // Update this node.
      clSize = 0;
      const uint16_t bitPosOfLastChunk = this->calcBitPos(2 * (childIdx + numChild_del + curNumAfterDel));
      const uint16_t bitSizeOfLastChunk = this->stccSize_ - bitPosOfLastChunk;
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
        this->reserveBitCapacity(sumWR);
        for (uint8_t i = 0; i < clSize; ++i) {
          this->stcc_.mvVals(this->getConstPtr_vals(), std::get<0>(changeList[i]), std::get<1>(changeList[i]), std::get<2>(changeList[i]));
        }
        this->stccSize_ = sumWR;
        this->shrinkBitCapacity();
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
    //   uint64_t sumWL = lnode->stccSize_;
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
    //     const auto w = this->stccSize_ - bitPos;
    //     changeList[clSize++] = {bitPos, sumWL, w};
    //     sumWL += w;
    //     lnode->stcc_.mvWCodes(this->getConstPtr_wCodes(), 2 * (childIdx + numChild_del), 2 * numL, 2 * numToLeft2);
    //   }
    //   lnode->numChildren_ = numL_new;
    //   lnode->updateWCodesAuxM(2 * numL_old, 2 * numL_new);
    //   { // Update vals of left node
    //     if (lnode->stccCapacity_ < sumWL) {
    //       lnode->stccCapacity_ = lnode->setBitCapacity(sumWL);
    //     }
    //     for (uint8_t i = 0; i < clSize; ++i) {
    //       lnode->mvVals(this->getConstPtr_vals(), std::get<0>(changeList[i]), std::get<1>(changeList[i]), std::get<2>(changeList[i]));
    //     }
    //     lnode->stccSize_ += sumWL;
    //   }

    //   if (numTailInR) {
        
    //   } else {
        
    //   }

    //   if (numSrcWCodesInL) {
    //     const uint64_t sumWL_ins = StepCodeUtil::sumW(srcWCodes, 0, 2 * numSrcWCodesInL);
    //     const uint64_t tailBitPos_new = this->calcBitPos(2 * childIdx) + sumWL_ins;
    //     this->stccSize_ = tailBitPos_new;
    //     const auto numTail = numL_new - (childIdx + numChild_ins);
    //     if (numTail) {
    //       const auto tailBitPos_old = this->calcBitPos(2 * (childIdx + numChild_del));
    //       const auto w = this->calcBitPos(2 * (childIdx + numChild_del + numTail)) - tailBitPos_old;
    //       this->stccSize_ += w;
    //       if (tailBitPos_new != tailBitPos_old) {
    //         if (this->stccCapacity_ < this->stccSize_) {
    //           this->stccCapacity_ = this->setBitCapacity(this->stccSize_);
    //         }
    //         this->mvVals(this->getConstPtr_vals(), tailBitPos_old, tailBitPos_new, w);
    //       }
    //       if (numChild_ins != numChild_del) {
    //         this->stcc_.mvWCodes(this->getConstPtr_wCodes(), 2 * (childIdx * numChild_del), 2 * (childIdx * numChild_ins), 2 * numTail);
    //       }
    //     } else {
    //       if (this->stccCapacity_ < this->stccSize_) {
    //         this->stccCapacity_ = this->setBitCapacity(this->stccSize_);
    //       }
    //     }
    //     this->stcc_.mvWCodes(srcWCodes, 0, 2 * childIdx, 2 * numSrcWCodesInL);
    //     this->updateWCodesAuxM(2 * childIdx, 2 * numL_new);
    //   } else {
    //     this->stccSize_ = this->calcBitPos(2 * numL_new); // shrink (just change bitSize)
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
      {//debug
        std::cerr << __func__ << ": numChild_ins = " << numChild_ins << ", childIdx = " << childIdx << ", numChild_del = " << numChild_del << std::endl;
      }
      assert(childIdx + numChild_del <= numChildren_);

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
        const uint16_t w = this->stccSize_ - bitPos;
        changeList[clSize++] = {bitPos, sumWR_increment, w};
        sumWR_increment += w;
        rnode->stcc_.mvWCodes(this->getConstPtr_wCodes(), 2 * (numL_old - numToRight2), 2 * numR_increment, 2 * numToRight2);
      }
      rnode->numChildren_ = static_cast<uint8_t>(numR_new);
      rnode->updateWCodesAuxM(0, 2 * numR_new);
      { // Update vals of "rnode".
        rnode->reserveBitCapacity(rnode->stccSize_ + sumWR_increment);
        if (numR_old) {
          rnode->stcc_.mvVals(rnode->getConstPtr_vals(), 0, sumWR_increment, rnode->stccSize_);
        }
        for (uint8_t i = 0; i < clSize; ++i) {
          rnode->stcc_.mvVals(this->getConstPtr_vals(), std::get<0>(changeList[i]), std::get<1>(changeList[i]), std::get<2>(changeList[i]));
        }
        rnode->stccSize_ += sumWR_increment;
      }

      if (numSrcWCodesInL) {
        const uint16_t sumWL_ins = static_cast<uint16_t>(StepCodeUtil::sumW(srcWCodes, 0, 2 * numSrcWCodesInL));
        const uint16_t tailBitPos_new = this->calcBitPos(2 * childIdx) + sumWL_ins;
        this->stccSize_ = tailBitPos_new;
        const uint16_t numTail = numL_new - (childIdx + numSrcWCodesInL);
        if (numTail) {
          const uint16_t tailBitPos_old = this->calcBitPos(2 * (childIdx + numChild_del));
          const uint16_t w = this->calcBitPos(2 * (childIdx + numChild_del + numTail)) - tailBitPos_old;
          this->stccSize_ += w;
          if (tailBitPos_new != tailBitPos_old) {
            this->reserveBitCapacity(this->stccSize_);
            this->stcc_.mvVals(this->getConstPtr_vals(), tailBitPos_old, tailBitPos_new, w);
          }
          this->stcc_.mvWCodes(this->getConstPtr_wCodes(), 2 * (childIdx + numChild_del), 2 * (childIdx + numChild_ins), 2 * numTail);
        } else {
          this->reserveBitCapacity(this->stccSize_);
        }
        this->stcc_.mvWCodes(srcWCodes, 0, 2 * childIdx, 2 * numSrcWCodesInL);
        this->updateWCodesAuxM(2 * childIdx, 2 * numL_new);
      } else { // shrink
        this->stccSize_ = this->calcBitPos(2 * numL_new); // shrink (just change bitSize)
        this->shrinkBitCapacity();
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
     const uint16_t numChild_ins,
     const uint16_t numChild_del //!< Length of wCodes of tgt to delete.
     ) noexcept {
      {//debug
        std::cerr << __func__ << ": childIdx = " << childIdx << ", sumW_ins = " << (int)sumW_ins << ", sumW_del = " << sumW_del
                  << ", numChild_ins = " << (int)numChild_ins << ", numChild_del = " << (int)numChild_del << std::endl;
      }
      const uint16_t tailNum = this->numChildren_ - (childIdx + numChild_del); // at least 0 by assumption.
      uint16_t tailW = 0;
      if (tailNum) {
        tailW = this->stccSize_ - this->calcBitPos(2 * (childIdx + numChild_del));
        this->stcc_.mvWCodes(this->stcc_.getConstPtr_wCodes(), 2 * (childIdx + numChild_del), 2 * (childIdx + numChild_ins), 2 * tailNum);
      }
      this->stcc_.mvWCodes(srcWCodes, 0, 2 * childIdx, 2 * numChild_ins);
      this->numChildren_ += numChild_ins - numChild_del;
      this->updateWCodesAuxM(2 * childIdx, 2 * this->numChildren_);
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
     const uint16_t numChild_ins
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
     uint16_t childIdx_ins, //!< Beginning childIdx relative from beginning position of "this" node
     const uint64_t * newWeights, //!< Storing weights to insert to stcc
     const uint64_t * newVals, //!< Storing new values to insert
     uint16_t numChild_ins
     ) {
      auto lnode = this;
      uint16_t numL = lnode->getNumChildren();
      uint16_t numNewElemInL = 0;
      if (childIdx_ins < numL) { // Insert new elements to lnode
        uint64_t bitPos = lnode->calcBitPos(2 * childIdx_ins);
        numNewElemInL = std::min(static_cast<uint16_t>(numChild_ins), static_cast<uint16_t>(numL - childIdx_ins));
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
        childIdx_ins += numNewElemInL - numL;
        uint64_t bitPos = rnode->calcBitPos(2 * childIdx_ins);
        for (uint16_t ci = childIdx_ins; ci < childIdx_ins + numChild_ins - numNewElemInL; ++ci) {
          uint8_t w = rnode->stcc_.readW(2 * ci);
          rnode->stcc_.writeWBits(newWeights[ci - childIdx_ins + numNewElemInL], bitPos, w);
          bitPos += w;
          w = rnode->stcc_.readW(2 * ci + 1);
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
     const uint16_t sumW_ins,
     const uint16_t sumW_del,
     const uint16_t numChild_ins,
     const uint16_t numChild_del //!< Length of wCodes of tgt to delete.
     ) noexcept {
      {//debug
        std::cerr << __func__ << ": childIdx = " << childIdx << ", sumW_ins = " << sumW_ins << ", sumW_del = " << sumW_del
                  << ", numChild_ins = " << (int)numChild_ins << ", numChild_del = " << (int)numChild_del
                  << ", this(" << this << ")->getNumChildren() = " << this->getNumChildren()
                  << ", rnode(" << rnode << ")->getNumChildren() = " << rnode->getNumChildren() << std::endl;
        // this->printStatistics(std::cout, true);
        // rnode->printStatistics(std::cout, true);
      }
      const uint16_t sumW_diff = static_cast<uint16_t>(sumW_ins - sumW_del);
      uint16_t growStccSize = this->stccSize_;
      uint16_t growNum = this->numChildren_;
      {
        const uint16_t newBitSize = static_cast<uint16_t>(growStccSize + rnode->stccSize_ + sumW_diff);
        this->reserveBitCapacity(newBitSize);
        this->stccSize_ = newBitSize;
      }
      if (childIdx < growNum) {
        const uint16_t tailNum = growNum - (childIdx + numChild_del); // at least 0 by assumption.
        if (tailNum) {
          const uint16_t tailW = growStccSize - this->calcBitPos(2 * (childIdx + numChild_del));
          this->stcc_.mvWCodes(this->stcc_.getConstPtr_wCodes(), 2 * (childIdx + numChild_del), 2 * (childIdx + numChild_ins), 2 * tailNum);
          this->stcc_.mvVals(this->stcc_.getConstPtr_vals(), growStccSize - tailW, static_cast<uint16_t>(growStccSize + sumW_diff - tailW), tailW);
        }
        this->stcc_.mvWCodes(srcWCodes, 0, 2 * childIdx, 2 * numChild_ins);
        growNum += numChild_ins - numChild_del;
        growStccSize += sumW_diff;
        this->stcc_.mvWCodes(rnode->stcc_.getConstPtr_wCodes(), 0, 2 * growNum, 2 * rnode->numChildren_);
        this->stcc_.mvVals(rnode->stcc_.getConstPtr_vals(), 0, growStccSize, rnode->stccSize_);
        this->numChildren_ = static_cast<uint8_t>(growNum + rnode->numChildren_);
      } else {
        const uint16_t childIdxInR = childIdx - growNum;
        if (childIdxInR) {
          const uint16_t stccSize1 = rnode->calcBitPos(2 * childIdxInR);
          this->stcc_.mvWCodes(rnode->stcc_.getConstPtr_wCodes(), 0, 2 * growNum, 2 * childIdxInR);
          this->stcc_.mvVals(rnode->stcc_.getConstPtr_vals(), 0, growStccSize, stccSize1);
          growNum += childIdxInR;
          growStccSize += stccSize1;
        }
        this->stcc_.mvWCodes(srcWCodes, 0, 2 * growNum, 2 * numChild_ins);
        growNum += numChild_ins;
        growStccSize += sumW_ins;
        const uint16_t tailNum = rnode->numChildren_ - (childIdxInR + numChild_del);
        if (tailNum) {
          const uint16_t tailW = rnode->stccSize_ - rnode->calcBitPos(2 * (childIdxInR + numChild_del));
          this->stcc_.mvWCodes(rnode->stcc_.getConstPtr_wCodes(), 2 * (childIdxInR + numChild_del), 2 * growNum, 2 * tailNum);
          this->stcc_.mvVals(rnode->stcc_.getConstPtr_vals(), rnode->stccSize_ - tailW, growStccSize, tailW);
        }
        this->numChildren_ = static_cast<uint8_t>(growNum + tailNum);
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
      {//debug
        std::cerr << __func__ << ": numChild_ins = " << numChild_ins << ", childIdx = " << childIdx << ", numChild_del = " << numChild_del << std::endl;
      }
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
        makeSpace(childIdx, wCodesTemp, sumW_ins, sumW_del, numChild_ins, numChild_del);
        this->writeNewElem(childIdx, weights, vals, numChild_ins);
        return;
      }

      if (idxInSibling) { // Check previous sibling.
        auto lnode = reinterpret_cast<BtmNodeT *>(parent->getChildPtr(idxInSibling - 1));
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
      {//debug
        std::cerr << __func__ << "(" << this << "): numChild_ins = " << numChild_ins
                  << ", childIdx = " << childIdx << ", numChild_del = " << numChild_del
                  << ", numChildren_ = " << (int)numChildren_
                  << ", stccSize _ = " << stccSize_ << ", stccCapacity_ = " << stccCapacity_
                  << ", parent = " << parent << ", idxInSibling = " << (int)idxInSibling << std::endl;
      }
      assert(numChild_ins <= numChild_del); // Assume that # insertion is at most # deletion.
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
        parent->changePSumFrom(idxInSibling, static_cast<int64_t>(sumWeights_ins - sumWeights_del));
      }

      const uint16_t num_new = numChildren_ - numChild_del + numChild_ins;
      if (num_new == 0) {
        delete this;
        parent->handleDeleteOfBtm(idxInSibling);
        return;
      }
      if (num_new < kBtmB / 2) {
        if (idxInSibling) { // Check previous sibling if it can be merged
          auto lnode = reinterpret_cast<BtmNodeT *>(parent->getChildPtr(idxInSibling - 1));
          const uint16_t lnum = lnode->getNumChildren();
          if (lnum + num_new <= 3 * kBtmB / 4) {
            lnode->merge(this, lnum + childIdx, wCodesTemp, sumW_ins, sumW_del, numChild_ins, numChild_del);
            lnode->updateWCodesAuxM(2 * lnum, 2 * lnode->numChildren_);
            if (numChild_ins) {
              lnode->writeNewElem(lnum + childIdx, weights, vals, numChild_ins);
            }
            parent->changePSumAt(idxInSibling - 1, parent->getPSum(idxInSibling + 1));
            delete this;
            parent->handleDeleteOfBtm(idxInSibling);
            return;
          }
        }
        if (idxInSibling + 1 < parent->getNumChildren()) { // Check next sibling if it can be merged
          auto rnode = reinterpret_cast<BtmNodeT *>(parent->getChildPtr(idxInSibling + 1));
          const uint16_t rnum = rnode->getNumChildren();
          if (rnum + num_new <= 3 * kBtmB / 4) {
            this->merge(rnode, childIdx, wCodesTemp, sumW_ins, sumW_del, numChild_ins, numChild_del);
            this->updateWCodesAuxM(2 * childIdx, 2 * this->numChildren_);
            if (numChild_ins) {
              this->writeNewElem(childIdx, weights, vals, numChild_ins);
            }
            parent->changePSumAt(idxInSibling, parent->getPSum(idxInSibling + 2));
            delete rnode;
            parent->handleDeleteOfBtm(idxInSibling + 1);
            return;
          }
        }
      }

      makeSpace(childIdx, wCodesTemp, sumW_ins, sumW_del, numChild_ins, numChild_del);
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
      {//debug
        std::cerr << __func__ << ": num = " << num << ", idx = " << idx << std::endl;
      }
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
        if (i % 2 == 0) {
          weights_diff += array[i - idx] - stcc_.readWBits(bitPos, w_old);
        }
        stcc_.writeWCode(StepCodeUtil::calcWCodeFromSteppedW(w_new), i);
        bitPos += w_old;
      }
      this->updateWCodesAuxM(idx, idx + num);

      if (sumW_ins != sumW_del) {
        const uint16_t newStccSize = this->stccSize_ + sumW_ins - sumW_del;
        if (newStccSize != this->stccSize_) {
          this->reserveBitCapacity(newStccSize);
          this->stcc_.mvVals(this->stcc_.getConstPtr_vals(), bitPos, bitPos0 + sumW_ins, this->stccSize_ - bitPos);
        }
        this->stccSize_ = newStccSize;
      }
      bitPos = bitPos0;
      for (uint16_t i = idx; i < idx + num; ++i) {
        uint8_t w = this->stcc_.readW(i);
        this->stcc_.writeWBits(array[i - idx], bitPos, w);
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
      return size + stccCapacity_ / 8;
    }


    size_t calcMemBytesDynArray() const noexcept {
      return stccCapacity_ / 8;
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
      os << "bit size = " << stccSize_ << ", bit capacity = " << stccCapacity_ << std::endl;
      os << "Total: " << calcMemBytes() << " bytes" << std::endl;
      os << "Memory usage for dynamic arrays of step code: " << calcMemBytesDynArray() << " bytes" << std::endl;
      if (verbose) {
        {
          const uint16_t num = static_cast<uint16_t>(numChildren_) * 2;
          os << "dump bit witdth stored in wCodes (" << stcc_.getConstPtr_wCodes() << ")" << std::endl;
          for (uint64_t i = 0; i < num; ++i) {
            os << (uint64_t)(stcc_.readW(i)) << " ";
          }
          os << std::endl;
        }
        {
          const uint16_t num = static_cast<uint16_t>(numChildren_) * 2;
          os << "dump values" << std::endl;
          for (uint64_t i = 0; i < num; ++i) {
            os << stcc_.read(i) << " ";
          }
          os << std::endl;
        }
        {
          os << "dump bits in vals_ (" << stcc_.getConstPtr_vals() << ")" << std::endl;
          for (uint64_t i = 0; i < (stccSize_ + 63) / 64; ++i) {
            os << "(" << i << ")";
            for (uint64_t j = 0; j < 64; ++j) {
              os << bits::readWBits_S(stcc_.getConstPtr_vals(), 64 * i + 63 - j, ctcbits::UINTW_MAX(1));
            }
            os << " ";
          }
          os << std::endl;
        }
      }
      os << "BtmNodeForPSumWithVal object (" << this << ") " << __func__ << "(" << verbose << ") END" << std::endl;
    }
  };




  /*!
   * @brief Dynamic partial sum data structure implemented by B+tree, where each leaf is associated with value.
   * @tparam kB Arity for internal node of B+tree, which should be in {32, 64, 128}.
   * @tparam kBtmB Arity for bottom node of B+tree, which should be in {16, 32, 64, 128}.
   * @par Notation
   */
  template<typename tparam_BTreeNodeT, typename tparam_BtmNodeT>
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
      return (sr_.root_ != nullptr);
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


    const BTreeNodeT * getConstRoot() const noexcept {
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
      assert(pos < getSumOfWeight());

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
     const bool verbose
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
        if (verbose) {
          getConstRoot()->printStatistics(os, true);
          os << "Btm Nodes" << std::endl;
          uint8_t idxInSib = 0;
          auto borderNode = sr_.root_->getLmBorderNode_DirectJump();
          while (reinterpret_cast<uintptr_t>(borderNode) != BTreeNodeT::NOTFOUND) {
            reinterpret_cast<BtmNodeT *>(borderNode->getChildPtr(idxInSib))->printStatistics(os, true);
            borderNode = borderNode->getNextBtmRef_DirectJump(idxInSib);
          }
        }
      } else {
        os << "Data structure is empty (not ready)." << std::endl;
      }
      os << "PSumWithValue object (" << this << ") " << __func__ << "(" << verbose << ") END" << std::endl;
    }


    uint64_t debugBtm
    (
     const BtmNodeT * btm,
     std::ostream & os
     ) const noexcept {
      const bool isLmBtm = (btm == reinterpret_cast<const BtmNodeT *>(getConstRoot()->getLmBtm_NoDirectJump()));
      const uint64_t sumWeights = btm->calcSumOfWeight();
      const uint16_t numChildren = btm->getNumChildren();
      uint64_t sum = 0;
      for (uint16_t i = isLmBtm; i < numChildren; ++i) {
        const uint64_t weight = btm->getWeight(i);
        if (weight == 0) {
          os << "error!! weightless element: i = " << (int)i << ", this = " << this << std::endl;
        }
        sum += weight;
      }
      if (sum != sumWeights) {
        os << "error!! weights does not match: sum = " << sum << ", sumWeights = " << sumWeights << std::endl;
      }
      return sumWeights;
    }


    void debugNode
    (
     const BTreeNodeT * node,
     std::ostream & os
     ) const noexcept {
      const uint16_t numChildren = node->getNumChildren();
      if (node->isBorder()) {
        for (uint16_t i = 0; i < numChildren; ++i) {
          const auto btm = reinterpret_cast<BtmNodeT *>(node->getChildPtr(static_cast<uint8_t>(i)));
          const uint64_t btmWeight = debugBtm(btm, os);
          const uint64_t psumWeight = node->getWeightOfChild(static_cast<uint8_t>(i));
          if (btmWeight != psumWeight) {
            os << "error!! btm weight: i = " << (int)i << ", btmWeight = " << btmWeight << ", psumWeight = " << psumWeight << std::endl;
            node->printStatistics(os, true);
            btm->printStatistics(os, true);
          }
        }
      } else {
        for (uint16_t i = 0; i < numChildren; ++i) {
          const auto child = node->getChildPtr(static_cast<uint8_t>(i));
          const uint64_t childWeight = child->getSumOfWeight();
          const uint64_t psumWeight = node->getWeightOfChild(static_cast<uint8_t>(i));
          if (childWeight != psumWeight) {
            os << "error!! weight: i = " << (int)i << ", childWeight = " << childWeight << ", psumWeight = " << psumWeight << std::endl;
            node->printStatistics(os, true);
            child->printStatistics(os, true);
          }
          if (node != child->getParent()) {
            os << "error!! child-parent: i = " << (int)i << ", parent = " << node << ", parent of child = " << child->getParent() << std::endl;
            node->printStatistics(os, true);
            child->printStatistics(os, true);
          }
          debugNode(child, os);
        }
      }
    }


    void printDebugInfo
    (
     std::ostream & os
     ) const noexcept {
      os << __func__ << ": PSumWithValue" << std::endl;
      if (!isReady()) {
        return;
      }
      // printStatistics(os, true);
      debugNode(sr_.root_, os);
    }
  };
} // namespace itmmti

#endif
