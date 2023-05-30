// (C)  Copyright Alexander (Lex) Kemper akemper@ncsu.edu 2022


#include<vector>
#include<cassert>
#include<iostream>

#include "cartan/Pauli.h"
#include "cartan/dla_generator.h"

#define MAX_LEVEL 40

#ifdef _OPENMP
    #include<omp.h>
#endif

// Add: pass in a function for the involution
std::tuple<algebra, algebra> get_algebra_by_commuting(
        algebra* p_algebra_g, 
        algebra* p_algebra_h,
        const std::optional<involutionFunc>& involution)
{	
    algebra& algebra_g = *p_algebra_g;
    algebra& algebra_h = *p_algebra_h;

    bool have_involution = false;
    bool use_h_algebra = false;

    // Did I get passed a pointer to an h algebra?
    if(p_algebra_h != nullptr)
    {
        use_h_algebra = true;
    }
    
    // Did I get passed an involution?
    if(involution)
    {
        std::cout << "Using an involution" << std::endl;
        have_involution = true;
    }


    // Copy the set into a vector
    std::vector<Pauli> algebra_g_keys;
    for(const auto& p : algebra_g)
    {
        algebra_g_keys.push_back(p);
    }

    // Prepare objects for algebra_k and algebra_m
    algebra algebra_k, algebra_m;

    // Initial algebra is just the Hamiltonian -- put that in m
    if(have_involution)
    for(const auto& term : algebra_g_keys)
    {
        assert(not parityY(term.code));
        algebra_m.insert(term);
    }




    // Do all commutators
    size_t endpoint = algebra_g_keys.size();
    size_t previous_endpoint = algebra_g_keys.size();
    int level = 1;
    for( ; level < MAX_LEVEL; level++)
    {
        size_t terms_added = 0;
        std::cout << "Currently have " << algebra_g.size() << " terms." << std::endl;

        previous_endpoint = endpoint;
        endpoint = algebra_g_keys.size();

#pragma omp parallel
{
#ifdef _OPENMP
        int my_th_id = omp_get_thread_num();
        //printf("Hello from thread %i\n", my_th_id);
#endif
        long terms_checked = 0, comms_calculated = 0, criticals=0;
        #pragma omp for schedule(static, 1)
        for(size_t ix=0; ix < endpoint; ix++)
        {
            const pauli_int Pauli1 = algebra_g_keys[ix].code;
            
            size_t inner_loop_start;
            if(ix < previous_endpoint and level > 1)
                inner_loop_start = previous_endpoint;
            else
                inner_loop_start = ix+1;

            for(size_t iy=inner_loop_start; iy < endpoint; iy++)
            {
                const pauli_int Pauli2 = algebra_g_keys[iy].code;

                comms_calculated++;

                //std::cout << "Checking " << Pauli(Pauli1) << " " << Pauli(Pauli2) << std::endl;

                pauli_int comm = symplectic_binary_comm(Pauli1, Pauli2);
                if(comm==0)
                    continue;


                // Add to the algebra if it's not in there already
                
                //std::cout << "Found one! " << Pauli(comm) << std::endl;
                terms_checked++;
                if(not algebra_g.count(comm))
                {
                    criticals++;
                    bool is_in_k = parityY(comm);
                    #pragma omp critical
                    {
                    algebra_g.insert(comm);
                    algebra_g_keys.push_back(comm);
                    terms_added++;
                    if(have_involution)
                    if(is_in_k)
                        algebra_k.insert(comm);
                    else
                        algebra_m.insert(comm);
                    }
                }
                    
            }
        }
        /*
        printf("Thread %d calculated %d commutators and checked %d in the set, fraction %f. Critical counts: %d: %f\n",
                my_th_id,comms_calculated,terms_checked,double(terms_checked)/comms_calculated,criticals,
                double(criticals)/comms_calculated);
                */

} // End parallel region

        // Find the terms in m that should also go in h
        if(use_h_algebra)
            for(const auto& m : algebra_m)
            {
                bool commute_with_all = true;
                for(const auto& h : algebra_h)
                {
                    if(symplectic_binary_comm(m.code, h.code)) // Do they not commute?
                    {
                        commute_with_all = false;
                        break;
                    }
                }
                if(commute_with_all)
                    algebra_h.insert(m);
            }

        std::cout << "Terms added: " << terms_added << std::endl;
        std::cout << "Algebra sizes: " << std::endl;
        std::cout << "|g| = " << algebra_g.size() << std::endl;;
        std::cout << "|k| = " << algebra_k.size() << std::endl;;
        std::cout << "|m| = " << algebra_m.size() << std::endl;;
        if(use_h_algebra)
        {
            std::cout << "|h| = " << algebra_h.size() << std::endl;;
            std::cout << "|k| + |h| = " << algebra_k.size() + algebra_h.size() << std::endl;
        }
        std::cout << "|k| + |m| = " << algebra_k.size() + algebra_m.size() << std::endl;

        if(terms_added == 0)// or (use_h_algebra and (algebra_k.size() + algebra_h.size() == algebra_m.size())))
            break;



        std::cout << "Hamiltonian (level " << level << ")" << std::endl;
        std::cout << "  " << terms_added << " terms" << std::endl;
        std::cout << '\n' << std::endl;

    } // End level loop

    std::cout << "Found " << algebra_g.size() << " terms." << std::endl;
    std::cout << "Found " << algebra_g_keys.size() << " terms in keys.\n" << std::endl;
#ifdef _OPENMP
    if(algebra_g.size() != algebra_g_keys.size())
    {
        std::cout << "Note: a potential mismatch between the number of terms in keys and algebra is possible." << std::endl;
        std::cout << "This is due to a race condition that is rare, but more importantly, does not affect anything. " << std::endl;
        std::cout << "\n" << std::endl;
    }
#endif

    if(level == MAX_LEVEL-1)
        std::cout << "Finished because of hard coded level limit." << std::endl;

    return {algebra_k, algebra_m};
}
