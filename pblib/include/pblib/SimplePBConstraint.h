#ifndef SIMPLEPBCONSTRAINT_H
#define SIMPLEPBCONSTRAINT_H

#include <assert.h>
#include "pbconstraint.h"
#include <vector>


enum PBTYPE {DONTCARE, AMO, AMK, PB};

class SimplePBConstraint : public PBLib::PBConstraint
{
protected:
    PBTYPE type;
    int64_t max_sum;
    int64_t max_weight;

public:
    SimplePBConstraint(int64_t max_sum, int64_t max_weight, PBTYPE type, std::vector<PBLib::WeightedLit>  const & literals, PBLib::Comparator comparator, int64_t less_eq, int64_t greater_eq);
    SimplePBConstraint(int64_t max_sum, int64_t max_weight, PBTYPE type, std::vector<PBLib::WeightedLit> const & literals, PBLib::Comparator comparator, int64_t bound);
    ~SimplePBConstraint() override = default;

    PBTYPE getType() const;
    int64_t getMaxWeight() const;
    int64_t getMaxSum() const;

    void print(bool std_err = false) const override;
    void printNoNL(bool std_err = false) const override;
};

#endif // SIMPLEPBCONSTRAINT_H
