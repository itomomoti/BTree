#include "PSumWithValue.hpp"
#include <gtest/gtest.h>
#include <stdint.h>
#include <algorithm>
#include <list>
#include <random>

namespace itmmti
{
  class PSumWithValue_Test : public ::testing::Test
  {
  };


  template<class BTreeNodeT, class BtmNodeT>
  bool checkPSum
  (
   const BTreeNodeT * node
   ) {
    if (node->isBorder()) {
      for (uint8_t i = 0; i < node->getNumChildren(); ++i) {
        const auto btm = reinterpret_cast<BtmNodeT *>(node->getChildPtr(i));
        const uint64_t weight = btm->calcSumOfWeight();
        if (node->getWeightOfChild(i) != weight) {
          std::cout << (int)i << "-th child (" << btm << ") of node (" << node << "); btm node says " << weight << " and node says " << node->getWeightOfChild(i);
          return false;
        }
      }
    } else {
      for (uint8_t i = 0; i < node->getNumChildren(); ++i) {
        const auto child = node->getChildPtr(i);
        if (checkPSum<BTreeNodeT, BtmNodeT>(child) == false) {
          return false;
        }
        const uint64_t weight = child->getSumOfWeight();
        if (node->getWeightOfChild(i) != weight) {
          std::cout << (int)i << "-th child (" << child << ") of node (" << node << "); child says " << weight << " and node says " << node->getWeightOfChild(i);
          return false;
        }
      }
    }
    return true;
  }


  template<class PSumWithValueT>
  void printAll(PSumWithValueT & obj) {
    using BTreeNodeT = typename PSumWithValueT::BTreeNodeT;
    using BtmNodeT = typename PSumWithValueT::BtmNodeT;

    if (obj.isReady()) {
      uint8_t idx = 0;
      uint64_t sum = 0;

      // std::cout << "root: " << std::endl;
      // obj.getRoot()->printStatistics(std::cout, false);
      BTreeNodeT * node = obj.getRoot()->getLmBorderNode_DirectJump();
      while (reinterpret_cast<uintptr_t>(node) != BTreeNodeT::NOTFOUND) {
        BtmNodeT * btm = reinterpret_cast<BtmNodeT *>(node->getChildPtr(idx));
        std::cout << sum << "(" << btm << ")|>" << " ";
        for (uint8_t i = 0; i < btm->getNumChildren(); ++i) {
          std::cout << "(" << btm->getWeight(i) << ", " << btm->getVal(i) << ") ";
          sum += btm->getWeight(i);
        }
        node = node->getNextBtmRef_DirectJump(idx);
        std::cout << std::endl;
      }
      std::cout << "<|" << sum << std::endl;
    }
  }


  template<class PSumWithValueT>
  bool checkBitSize(PSumWithValueT & obj) {
    using BTreeNodeT = typename PSumWithValueT::BTreeNodeT;
    using BtmNodeT = typename PSumWithValueT::BtmNodeT;

    if (obj.isReady()) {
      uint8_t idx = 0;

      // std::cout << "root: " << std::endl;
      // obj.getRoot()->printStatistics(std::cout, false);
      BTreeNodeT * node = obj.getRoot()->getLmBorderNode_DirectJump();
      while (reinterpret_cast<uintptr_t>(node) != BTreeNodeT::NOTFOUND) {
        BtmNodeT * btm = reinterpret_cast<BtmNodeT *>(node->getChildPtr(idx));
        uint64_t sum = 0;
        for (uint8_t i = 0; i < btm->getNumChildren(); ++i) {
          sum += StepCodeUtil::calcSteppedW(btm->getWeight(i)) + StepCodeUtil::calcSteppedW(btm->getVal(i));
        }
        if (btm->getBitSize() != sum) {
          std::cout << "error: " << btm << ": " << "bitSize = " << btm->getBitSize() << ", sum = " << sum << std::endl;
          btm->printStatistics(std::cout, true);
          return false;
        }
        node = node->getNextBtmRef_DirectJump(idx);
      }
    }
    return true;
  }


  TEST_F(PSumWithValue_Test, Simple)
  {
    const size_t num = 32;

    const uint8_t kB = 32;
    const uint8_t kBtmB = 32;
    using BTreeNodeT = BTreeNode<kB>;
    using BtmNodeT = BtmNodeForPSumWithVal<kBtmB>;

    PSumWithValue<kB, kBtmB> obj;
    obj.init();
    obj.printStatistics(std::cout);

    {
      BTreeNodeT * node = nullptr;
      uint8_t idx;
      BtmNodeT * btmNode = obj.getRmBtm(node, idx);
      node->printStatistics(std::cout);
      btmNode->printStatistics(std::cout, true);
      {
        uint64_t weights[]  = {1, 2, 3, 4};
        uint64_t vals[] = {bits::UINTW_MAX(4), bits::UINTW_MAX(5), bits::UINTW_MAX(6), bits::UINTW_MAX(7)};
        // uint64_t weights[]  = {1, 2};
        // uint64_t vals[] = {bits::UINTW_MAX(2), bits::UINTW_MAX(4)};
        btmNode->insert(weights, vals, 4, 1, 0, node, idx);
        btmNode->insert(weights, vals, 4, 1, 0, node, idx);
        btmNode->insert(weights, vals, 4, 1, 0, node, idx);
        btmNode->insert(weights, vals, 4, 1, 0, node, idx);
        btmNode->insert(weights, vals, 4, 1, 0, node, idx);
        btmNode->insert(weights, vals, 4, 1, 0, node, idx);
        btmNode->insert(weights, vals, 4, 1, 0, node, idx);

        obj.printStatistics(std::cout);
        btmNode->printStatistics(std::cout, true);

        btmNode->insert(weights, vals, 4, 1, 0, node, idx);
        obj.printStatistics(std::cout);
        btmNode->printStatistics(std::cout, true);
      }
      node = node->getNextBtmRef_DirectJump(idx);
      BtmNodeT * btmNode2 = reinterpret_cast<BtmNodeT *>(node->getChildPtr(idx));
      btmNode2->printStatistics(std::cout, true);
      {
        uint64_t weights[]  = {5, 6, 7, 8, 9, 10};
        uint64_t vals[] = {bits::UINTW_MAX(9), bits::UINTW_MAX(10), bits::UINTW_MAX(9), bits::UINTW_MAX(10), bits::UINTW_MAX(9), bits::UINTW_MAX(10)};

        btmNode->insert(weights, vals, 6, 2, 0, node, idx-1);
        btmNode->printStatistics(std::cout, true);

        btmNode->insert(weights, vals, 6, 2, 0, node, idx-1);
        btmNode->printStatistics(std::cout, true);

        btmNode->insert(weights, vals, 6, 24, 0, node, idx-1);
        btmNode->printStatistics(std::cout, true);
        btmNode2->printStatistics(std::cout, true);

        btmNode->insert(weights, vals, 6, 24, 0, node, idx-1);
        btmNode->printStatistics(std::cout, true);
        btmNode2->printStatistics(std::cout, true);

        btmNode->insert(weights, vals, 6, 24, 0, node, idx-1);
        btmNode->printStatistics(std::cout, true);
        btmNode2->printStatistics(std::cout, true);

        btmNode->insert(weights, vals, 6, 24, 0, node, idx-1);
        btmNode->printStatistics(std::cout, true);
        btmNode2->printStatistics(std::cout, true);
      }
    }
  }


  TEST_F(PSumWithValue_Test, SimpleLarger)
  {
    const uint8_t kB = 32;
    const uint8_t kBtmB = 32;
    using BTreeNodeT = BTreeNode<kB>;
    using BtmNodeT = BtmNodeForPSumWithVal<kBtmB>;

    PSumWithValue<kB, kBtmB> obj;
    obj.init();

    uint64_t weights[128];
    uint64_t vals[128];
    for (uint8_t i = 0; i < 128; ++i) {
      weights[i] = 1;
      vals[i] = UINT64_C(1) << ((i % 16) * 4);
    }

    BTreeNodeT * node = nullptr;
    uint8_t idxInSibOfBtm;
    BtmNodeT * btm = obj.getRmBtm(node, idxInSibOfBtm);
    btm->insert(weights, vals, kBtmB - 1, 1, 0, node, idxInSibOfBtm);
    // node->printStatistics(std::cout, true);
    ++(weights[0]);
    btm->insert(weights, vals, kBtmB, kBtmB, 0, node, idxInSibOfBtm);
    // node->printStatistics(std::cout, true);

    {
      btm = reinterpret_cast<BtmNodeT *>(node->getChildPtr(1));
      size_t num = 8;
      // size_t num = 61;
      for (uint64_t j = 0; j < num; ++j) {
        ++(weights[0]);
        btm->insert(weights, vals, kBtmB, 0, 0, node, 1);
        // node->printStatistics(std::cout, true);
      }
      // uint64_t pos = obj.getSumOfWeight();
      // BTreeNodeT * node = nullptr;
      // uint8_t idxInSibOfBtm;
      // searchBtm(obj.getSumOfWeight(), node, idxInSibOfBtm);
    }
    obj.getRoot()->printStatistics(std::cout, false);
    obj.getRoot()->getChildPtr(0)->printStatistics(std::cout, false);
    obj.getRoot()->getChildPtr(1)->printStatistics(std::cout, false);
    printAll(obj);
  }


  TEST_F(PSumWithValue_Test, RandomInsert)
  {
    const uint8_t kB = 8;
    const uint8_t kBtmB = 8;
    using BTreeNodeT = BTreeNode<kB>;
    using BtmNodeT = BtmNodeForPSumWithVal<kBtmB>;

    std::list<uint64_t> naive_list;

    PSumWithValue<kB, kBtmB> obj;
    obj.init();
    obj.printStatistics(std::cout);

    uint64_t weights[128];
    uint64_t vals[128];
    for (uint8_t i = 0; i < 128; ++i) {
      weights[i] = 1;
      vals[i] = UINT64_C(1) << ((i % 16) * 4);
    }

    // std::random_device rnd; // 非決定的な乱数生成器を生成
    // std::mt19937 mt(rnd()); // メルセンヌ・ツイスタの32ビット版、引数は初期シード値
    std::mt19937 mt(0); // メルセンヌ・ツイスタの32ビット版、引数は初期シード値

    {
      // size_t num = 98;
      size_t num = 50000;
      for (uint64_t j = 0; j < num; ++j) {
        std::cout << "insert_loop: " << j << std::endl;

        {
          const uint64_t sum = obj.getSumOfWeight();
          uint64_t pos = 0;
          if (sum) {
            pos = (mt() % sum);
          }
          std::cout << "pos = " << pos << ", sum = " << sum;
          // if (j == 94) {
          //   std::cout << "before" << std::endl;
          //   BTreeNodeT * dummynode;
          //   uint8_t dummyidx;
          //   obj.getLmBtm(dummynode, dummyidx)->printStatistics(std::cout, true);
          // }
          auto itr = naive_list.begin();
          for (uint64_t i = 0; i < pos; ++i) {
            ++itr;
          }
          BTreeNodeT * node = nullptr;
          uint8_t idxInSibOfBtm;
          BtmNodeT * btm;
          uint8_t idx;
          if (pos < sum) {
            btm = obj.searchBtm(pos, node, idxInSibOfBtm);
            idx = btm->searchPos(pos);
          } else {
            btm = obj.getRmBtm(node, idx);
            ++idx;
          }
          const uint16_t num_ins = (mt() % kB) + 1;
          const uint16_t num_del = mt() % std::min((btm->getNumChildren() - idx + 1), 3);

          std::cout << ", btm = " << btm << ", num_ins = " << num_ins << ", idx = " << (int)idx << ", num_del = " << num_del << ", idxInSibOfBtm = " << (int)idxInSibOfBtm << std::endl;

          for (uint16_t i = 0; i < num_del; ++i) {
            itr = naive_list.erase(itr);
          }
          for (uint16_t i = 0; i < num_ins; ++i) {
            naive_list.insert(itr, vals[i]);
          }

          if (j == 599) {
            std::cout << "before" << std::endl;
            node->printStatistics(std::cout, true);
            btm->printStatistics(std::cout, true);
          }
          btm->insert(weights, vals, num_ins, idx, num_del, node, idxInSibOfBtm);
          if (j == 599) {
            std::cout << "after" << std::endl;
            node->printStatistics(std::cout, true);
            btm->printStatistics(std::cout, true);
            obj.getRoot()->printStatistics(std::cout, true);
          }
        }

        { // check correctness
          if (checkBitSize(obj) == false) {
            abort();
          }
          // std::cout << "checkBitSize passed" << std::endl;
          if (checkPSum<BTreeNodeT, BtmNodeT>(obj.getRoot()) == false) {
            abort();
          }
          // std::cout << "checkPSum passed" << std::endl;
          BTreeNodeT * node = nullptr;
          uint8_t idxInSibOfBtm;
          BtmNodeT * btm = obj.getLmBtm(node, idxInSibOfBtm);
          uint16_t idx = 0;
          uint64_t num = 0;
          for (auto itr = naive_list.begin(); itr != naive_list.end(); ++itr) {
            if (++idx >= btm->getNumChildren()) {
              node = node->getNextBtmRef_DirectJump(idxInSibOfBtm);
              btm = reinterpret_cast<BtmNodeT *>(node->getChildPtr(idxInSibOfBtm));
              idx = 0;
            }
            ++num;
            // if (252 <= j && j <= 254) {
            //   std::cout << "[" << num << "]" << *itr << ", ";
            // }
            ASSERT_EQ(*itr, btm->getVal(idx));
          }
          if (598 <= j && j <= 599) {
            std::cout << "<|" << num << std::endl;
          }
        }

        if (598 <= j && j <= 599) {
          obj.getRoot()->printStatistics(std::cout, true);
          printAll(obj);
        }
      }
      obj.getRoot()->printStatistics(std::cout, true);
      obj.printStatistics(std::cout, true);
    }
  }



} // namespace itmmti



