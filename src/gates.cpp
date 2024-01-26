#include "Simulator.h"
#include "util_sim.h"


void Simulator::Toffoli(int targ, std::vector<int> cont, std::vector<int> ncont, int cur_bit)
{
    assert((cont.size() + ncont.size()) < n);
    int IsBadtarg = 0;
    int cont_tot = cont.size() + ncont.size();
    for (int i = 0; i < cont_tot; i++)
    {
        if (i < cont.size())
        {
            if (targ == cont[i])
            {
                IsBadtarg = 1;
                break;
            }
        }
        else
        {
            if (targ == ncont[i - cont.size()])
            {
                IsBadtarg = 1;
                break;
            }
        }
    }
    assert(!IsBadtarg);
    DdNode *term1, *term2, *term3, *g, *tmp;

    g = Cudd_ReadOne(manager);
    Cudd_Ref(g);
    for (int h = cont.size() - 1; h >= 0; h--)
    {
        tmp = Cudd_bddAnd(manager, Cudd_bddIthVar(manager, cont[h]), g);
        Cudd_Ref(tmp);
        Cudd_RecursiveDeref(manager, g);
        g = tmp;
    }
    for (int h = ncont.size() - 1; h >= 0; h--)
    {
        tmp = Cudd_bddAnd(manager, Cudd_Not(Cudd_bddIthVar(manager, ncont[h])), g);
        Cudd_Ref(tmp);
        Cudd_RecursiveDeref(manager, g);
        g = tmp;
    }

            //term1
            term1 = Cudd_ReadOne(manager);
            Cudd_Ref(term1);
            tmp = Cudd_bddAnd(manager, Single_Bdd[cur_bit], term1);
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(manager, term1);
            term1 = tmp;
            tmp = Cudd_bddAnd(manager, Cudd_Not(g), term1);
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(manager, term1);
            term1 = tmp;

            //term2
            term2 = Cudd_Cofactor(manager, Single_Bdd[cur_bit], Cudd_Not(Cudd_bddIthVar(manager, targ)));
            Cudd_Ref(term2);

            tmp = Cudd_Cofactor(manager, term2, g);
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(manager, term2);
            term2 = tmp;

            tmp = Cudd_bddAnd(manager, term2, g);
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(manager, term2);
            term2 = tmp;

            tmp = Cudd_bddAnd(manager, term2, Cudd_bddIthVar(manager, targ));
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(manager, term2);
            term2 = tmp;

            //term3
            term3 = Cudd_Cofactor(manager, Single_Bdd[cur_bit], Cudd_bddIthVar(manager, targ));
            Cudd_Ref(term3);
            Cudd_RecursiveDeref(manager, Single_Bdd[cur_bit]);

            tmp = Cudd_Cofactor(manager, term3, g);
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(manager, term3);
            term3 = tmp;

            tmp = Cudd_bddAnd(manager, term3, g);
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(manager, term3);
            term3 = tmp;

            tmp = Cudd_bddAnd(manager, term3, Cudd_Not(Cudd_bddIthVar(manager, targ)));
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(manager, term3);
            term3 = tmp;

            //OR
            Single_Bdd[cur_bit] = Cudd_bddOr(manager, term1, term2);
            Cudd_Ref(Single_Bdd[cur_bit]);
            Cudd_RecursiveDeref(manager, term1);
            Cudd_RecursiveDeref(manager, term2);
            tmp = Cudd_bddOr(manager, term3, Single_Bdd[cur_bit]);
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(manager, term3);
            Cudd_RecursiveDeref(manager, Single_Bdd[cur_bit]);
            Single_Bdd[cur_bit] = tmp;
    Cudd_RecursiveDeref(manager, g);
    gatecount++;
    nodecount();
}

void Simulator::Fredkin(int swapA , int swapB, std::vector<int> cont, int cur_bit)
{
    assert(cont.size() < n);
    int IsBadtarg = 0;
    for (int i = 0; i < cont.size(); i++)
    {
        if ((swapA == cont[i]) || (swapB == cont[i]))
        {
            IsBadtarg = 1;
            break;
        }
    }
    assert(!IsBadtarg);

    DdNode *term1, *term2, *term3, *g, *tmp, *tmp0;

    g = Cudd_ReadOne(manager);
    Cudd_Ref(g);
    for (int h = cont.size() - 1; h >= 0; h--)
    {
        tmp = Cudd_bddAnd(manager, Cudd_bddIthVar(manager, cont[h]), g);
        Cudd_Ref(tmp);
        Cudd_RecursiveDeref(manager, g);
        g = tmp;
    }

            //term1
            term1 = Cudd_ReadOne(manager);
            Cudd_Ref(term1);
            tmp = Cudd_bddAnd(manager, Single_Bdd[cur_bit], term1);
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(manager, term1);
            term1 = tmp;
            tmp = Cudd_bddXor(manager, Cudd_bddIthVar(manager, swapA), Cudd_bddIthVar(manager, swapB));
            Cudd_Ref(tmp);
            tmp0 = Cudd_Not(Cudd_bddAnd(manager, g, tmp));
            Cudd_Ref(tmp0);
            Cudd_RecursiveDeref(manager, tmp);
            tmp = Cudd_bddAnd(manager, term1, tmp0);
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(manager, term1);
            Cudd_RecursiveDeref(manager, tmp0);
            term1 = tmp;

            //term2
            term2 = Cudd_Cofactor(manager, Single_Bdd[cur_bit], g);
            Cudd_Ref(term2);

            tmp = Cudd_Cofactor(manager, term2, Cudd_Not(Cudd_bddIthVar(manager, swapA)));
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(manager, term2);
            term2 = tmp;

            tmp = Cudd_Cofactor(manager, term2, Cudd_bddIthVar(manager, swapB));
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(manager, term2);
            term2 = tmp;

            tmp = Cudd_bddAnd(manager, term2, g);
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(manager, term2);
            term2 = tmp;

            tmp = Cudd_bddAnd(manager, term2, Cudd_bddIthVar(manager, swapA));
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(manager, term2);
            term2 = tmp;

            tmp = Cudd_bddAnd(manager, term2, Cudd_Not(Cudd_bddIthVar(manager, swapB)));
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(manager, term2);
            term2 = tmp;

            //term3
            term3 = Cudd_Cofactor(manager, Single_Bdd[cur_bit], g);
            Cudd_Ref(term3);

            tmp = Cudd_Cofactor(manager, term3, Cudd_bddIthVar(manager, swapA));
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(manager, term3);
            term3 = tmp;

            tmp = Cudd_Cofactor(manager, term3, Cudd_Not(Cudd_bddIthVar(manager, swapB)));
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(manager, term3);
            term3 = tmp;

            tmp = Cudd_bddAnd(manager, term3, g);
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(manager, term3);
            term3 = tmp;

            tmp = Cudd_bddAnd(manager, term3, Cudd_Not(Cudd_bddIthVar(manager, swapA)));
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(manager, term3);
            term3 = tmp;

            tmp = Cudd_bddAnd(manager, term3, Cudd_bddIthVar(manager, swapB));
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(manager, term3);
            term3 = tmp;

            //OR
            Single_Bdd[cur_bit] = Cudd_bddOr(manager, term1, term2);
            Cudd_Ref(Single_Bdd[cur_bit]);
            Cudd_RecursiveDeref(manager, term1);
            Cudd_RecursiveDeref(manager, term2);
            tmp = Cudd_bddOr(manager, term3, Single_Bdd[cur_bit]);
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(manager, term3);
            Cudd_RecursiveDeref(manager, Single_Bdd[cur_bit]);
            Single_Bdd[cur_bit] = tmp;
    Cudd_RecursiveDeref(manager, g);
    gatecount++;
    nodecount();
}

void Simulator::PauliX(int iqubit, int cur_bit)
{
    assert((iqubit >= 0) & (iqubit < n));

    DdNode *tmp, *term1, *term2;

            /*
            tmp = Cudd_bddCompose(manager,  Single_Bdd[cur_bit], Cudd_Not(Cudd_bddIthVar(manager, iqubit)), iqubit);
            Cudd_Ref(tmp);            
            
            Cudd_RecursiveDeref(manager, Single_Bdd[cur_bit]);
            Single_Bdd[cur_bit] = tmp;*/
            
            //term1
            term1 = Cudd_Cofactor(manager, Single_Bdd[cur_bit], Cudd_Not(Cudd_bddIthVar(manager, iqubit)));
            Cudd_Ref(term1);

            tmp = Cudd_bddAnd(manager, term1, Cudd_bddIthVar(manager, iqubit));
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(manager, term1);
            term1 = tmp;

            //term2
            term2 = Cudd_Cofactor(manager, Single_Bdd[cur_bit], Cudd_bddIthVar(manager, iqubit));
            Cudd_Ref(term2);
            Cudd_RecursiveDeref(manager, Single_Bdd[cur_bit]);

            tmp = Cudd_bddAnd(manager, term2, Cudd_Not(Cudd_bddIthVar(manager, iqubit)));
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(manager, term2);
            term2 = tmp;

            //OR
            Single_Bdd[cur_bit] = Cudd_bddOr(manager, term1, term2);
            Cudd_Ref(Single_Bdd[cur_bit]);
            Cudd_RecursiveDeref(manager, term1);
            Cudd_RecursiveDeref(manager, term2);
    gatecount++;
    nodecount();
}
