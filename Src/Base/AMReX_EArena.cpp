#ifdef CUDA

#include <utility>
#include <cstring>

#include <AMReX_EArena.H>

#include <AMReX_Device.H>
#include <AMReX_CUDA_helper.H>

#include <cuda_runtime_api.h>
#include <cuda.h>

namespace amrex {

EArena::EArena (size_t hunk_size)
{
    //
    // Force alignment of hunksize.
    //
    m_hunk = Arena::align(hunk_size == 0 ? DefaultHunkSize : hunk_size);
    m_used = 0;

    BL_ASSERT(m_hunk >= hunk_size);
    BL_ASSERT(m_hunk%Arena::align_size == 0);
}

EArena::~EArena ()
{
    CUresult cudaResult;
    CUcontext cuctx;
    cudaResult = cuCtxGetCurrent(&cuctx); 
    if (cudaResult != CUDA_ERROR_DEINITIALIZED) {
        for (unsigned int i = 0, N = m_alloc.size(); i < N; i++)
            checkCudaErrors(cudaFreeHost(m_alloc[i]));
    }
    // otherwise device has been reset in amrex::Finalize() and 
    // all device memory has been free. So we should not call cudaFree
}

void*
EArena::alloc_pinned (size_t nbytes)
{
    nbytes = EArena::align(nbytes == 0 ? 1 : nbytes);
    //
    // Find node in freelist at lowest memory address that'll satisfy request.
    //
    NL::iterator free_it = m_freelist.begin();

    for ( ; free_it != m_freelist.end(); ++free_it)
        if ((*free_it).size() >= nbytes)
            break;

    void* vp = 0;

    if (free_it == m_freelist.end())
    {
        const size_t N = nbytes < m_hunk ? m_hunk : nbytes;

        if ( (m_used + N) > max_heap_size )
            amrex::Abort("Not enough pinned memory available in EArena.");

        checkCudaErrors(cudaHostAlloc(&vp, N, cudaHostAllocDefault));

        m_used += N;

        m_alloc.push_back(vp);

        if (nbytes < m_hunk)
        {
            //
            // Add leftover chunk to free list.
            //
            // Insert with a hint -- should be largest block in the set.
            //
            void* block = static_cast<char*>(vp) + nbytes;

            m_freelist.insert(m_freelist.end(), Node(block, m_hunk-nbytes));
        }
    }
    else
    {
        BL_ASSERT((*free_it).size() >= nbytes);
        BL_ASSERT(m_busylist.find(*free_it) == m_busylist.end());

        vp = (*free_it).block();

        if ((*free_it).size() > nbytes)
        {
            //
            // Insert remainder of free block back into freelist.
            //
            // Insert with a hint -- right after the current block being split.
            //
            Node freeblock = *free_it;

            freeblock.size(freeblock.size() - nbytes);

            freeblock.block(static_cast<char*>(vp) + nbytes);

            m_freelist.insert(free_it, freeblock);
        }

        m_freelist.erase(free_it);
    }

    m_busylist.insert(Node(vp, nbytes));

    BL_ASSERT(!(vp == 0));

    return vp;
}

void
EArena::free_pinned (void* vp)
{
    if (vp == 0)
        //
        // Allow calls with NULL as allowed by C++ delete.
        //
        return;
    //
    // `vp' had better be in the busy list.
    //
    NL::iterator busy_it = m_busylist.find(Node(vp,0));

    BL_ASSERT(!(busy_it == m_busylist.end()));
    BL_ASSERT(m_freelist.find(*busy_it) == m_freelist.end());
    //
    // Put free'd block on free list and save iterator to insert()ed position.
    //
    std::pair<NL::iterator,bool> pair_it = m_freelist.insert(*busy_it);

    BL_ASSERT(pair_it.second == true);

    NL::iterator free_it = pair_it.first;

    BL_ASSERT(free_it != m_freelist.end() && (*free_it).block() == (*busy_it).block());
    //
    // And remove from busy list.
    //
    m_busylist.erase(busy_it);
    //
    // Coalesce freeblock(s) on lo and hi side of this block.
    //
    if (!(free_it == m_freelist.begin()))
    {
        NL::iterator lo_it = free_it;

        --lo_it;

        void* addr = static_cast<char*>((*lo_it).block()) + (*lo_it).size();

        if (addr == (*free_it).block())
        {
            bool merge = true;
            for (unsigned int i = 0, N = m_alloc.size(); i < N; i++) {
                // Don't merge two nodes if the merge will give a node 
                // whose memory block crosses the hunk boundary
                if (addr == m_alloc[i]) {
                    merge = false;
                }
            }
            //
            // This cast is needed as iterators to set return const values;
            // i.e. we can't legally change an element of a set.
            // In this case I want to change the size() of a block
            // in the freelist.  Since size() is not used in the ordering
            // relations in the set, this won't effect the order;
            // i.e. it won't muck up the ordering of elements in the set.
            // I don't want to have to remove the element from the set and
            // then reinsert it with a different size() as it'll just go
            // back into the same place in the set.
            //
            if (merge) {
                Node* node = const_cast<Node*>(&(*lo_it));
                BL_ASSERT(!(node == 0));
                node->size((*lo_it).size() + (*free_it).size());
                m_freelist.erase(free_it);
                free_it = lo_it;
            }
        }
    }

    NL::iterator hi_it = free_it;

    void* addr = static_cast<char*>((*free_it).block()) + (*free_it).size();

    if (++hi_it != m_freelist.end() && addr == (*hi_it).block())
    {
        bool merge = true;
        for (unsigned int i = 0, N = m_alloc.size(); i < N; i++) {
            // Don't merge two nodes if the merge will give a node 
            // whose memory block crosses the hunk boundary
            if (addr == m_alloc[i]) {
                merge = false;
            }
        }
        if (merge) {
            //
            // Ditto the above comment.
            //
            Node* node = const_cast<Node*>(&(*free_it));
            BL_ASSERT(!(node == 0));
            node->size((*free_it).size() + (*hi_it).size());
            m_freelist.erase(hi_it);
        }
    }
}

std::size_t
amrex::EArena::align (std::size_t s)
{
    std::size_t x = s + (align_size-1);
    x -= x & (align_size-1);
    return x;
}


size_t
EArena::heap_space_used () const
{
    return m_used;
}

}

#endif