#pragma once

#include <vfs/containers/EmptyArrayImpl.h>
#include <vfs/containers/VecImpl.h>
#include <vfs/support/Opt.h>

///
template<class Ci_nb_dims, class Ct_sval, class M_dirs, class M_offs, class B_dirs, class B_offs>
struct LegendreTransform {
    static constexpr int nb_dims              = Ci_nb_dims::value;
    using                Sval                 = Ct_sval::value;
    using                Svec                 = Vfs::VecImpl<Sval,nb_dims>;
    using                Ret                  = std::tuple<Vfs::VecImpl<Svec>,Vfs::VecImpl<Sval>,Vfs::VecImpl<Svec>,Vfs::VecImpl<Sval>>;
    using                PI                   = Vfs::PI;

    Ret                  transform_without_dir( Svec pos, Svec dir, bool add_bnd );
    Ret                  transform            ();

    void                 update_used_ms_and_bs( const auto &vertex_coords, const auto &vertex_cuts );
    void                 make_new_affs        ( Vfs::VecImpl<Svec> &nm_dirs, Vfs::VecImpl<Sval> &nm_offs, const auto &vertex_coords, const auto &vertex_cuts );
    void                 make_new_bnds        ( Vfs::VecImpl<Svec> &nb_dirs, Vfs::VecImpl<Sval> &nb_offs, const auto &vertex_coords, const auto &vertex_cuts );
    auto                 first_eq_bnd         () -> Vfs::Opt<std::pair<Svec,Svec>>;
    PI                   nb_inf_cuts          ( const auto *cuts );
    auto                 vertex_corr          ();
    auto                 unused_dir           () -> Vfs::Opt<std::pair<Svec,Svec>>;
    Vfs::Opt<PI>         aff_cut              ( PI cut ) const;
    Vfs::Opt<PI>         bnd_cut              ( PI cut ) const;
    Vfs::Opt<PI>         inf_cut              ( PI cut ) const;


    Ci_nb_dims           ci_nb_dims;
    Ct_sval              ct_sval;
    const M_dirs&        m_dirs;
    const M_offs&        m_offs;
    const B_dirs&        b_dirs;
    const B_offs&        b_offs;

    std::vector<bool>    used_ms;
    std::vector<bool>    used_bs;
};

/// return ( new_m_dirs, new_m_offs, new_b_dirs, new_b_offs )
auto legendre_transform( const auto &m_dirs, const auto &m_offs, const auto &b_dirs, const auto &b_offs, auto ci_nb_dims ) {
    // make the legendre transform
    auto ct_sval = DECAYED_CT_OF( m_dirs[ 0 ][ 0 ] + m_offs[ 0 ] + b_dirs[ 0 ][ 0 ] + b_offs[ 0 ] );
    LegendreTransform lt( ci_nb_dims, ct_sval, m_dirs, m_offs, b_dirs, b_offs );
    auto res = lt.transform();

    // (temporary) conversion to simple vectors. TODO: use list[list] => tensors
    Vfs::VecImpl<typename GET_DT_VALUE( ct_sval )> rm_dirs;
    for( const auto &v : std::get<0>( res ) )
        rm_dirs.append( v );

    Vfs::VecImpl<typename GET_DT_VALUE( ct_sval )> rb_dirs;
    for( const auto &v : std::get<2>( res ) )
        rb_dirs.append( v );

    return std::tuple( std::move( rm_dirs ), std::move( std::get<1>( res ) ), std::move( rb_dirs ), std::move( std::get<3>( res ) ) );
}

#include "legendre_transform.tcc" // IWYU pragma: export
