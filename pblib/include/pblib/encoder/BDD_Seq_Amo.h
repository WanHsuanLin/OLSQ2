/*
 *
 */

#ifndef BDD_SEQ_AMO_H
#define BDD_SEQ_AMO_H


#include <vector>
#include <sstream>

#include "../SimplePBConstraint.h"
#include "../IncSimplePBConstraint.h"
#include "../PBConfig.h"
#include "../clausedatabase.h"
#include "../auxvarmanager.h"
#include "../weightedlit.h"
#include "Encoder.h"

class BDD_Seq_Amo : public Encoder
{
private:
    std::vector<Lit> _literals;
    std::vector<Lit> aux;

public:
    void encode_intern( std::vector<Lit>& literals, ClauseDatabase & formula, AuxVarManager & auxvars);
    void encode(const SimplePBConstraint& pbconstraint, ClauseDatabase & formula, AuxVarManager & auxvars) override;
    int64_t encodingValue(const SimplePBConstraint& pbconstraint) override;

    BDD_Seq_Amo(PBConfig & config);
    ~BDD_Seq_Amo() override = default;
};

#endif // BDD_SEQ_AMO_H
