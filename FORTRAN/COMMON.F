      COMMON /elmatx/ x(3), y(3)
      COMMON /matl/ dxe, qe
      COMMON /bcs/ nbc, ib(maxbc), bv(maxbc)
      COMMON /av/ a(maxav), jgf, jgsm, np, nbw
      COMMON /nubsq/ erel(maxne), ebsq(maxne)
      COMMON /bh/ b(maxbh), h(maxbh)
      COMMON /newt/ esm(3,3), ef(maxne,3), esmuni(maxne,3,3)
      COMMON /nldat/ rel(maxbh), bsq(maxbh), const(3), grdmax
      COMMON /cgrad/ ns(3), nel(maxne,3), nmtl(maxne),
     + xc(maxnp), yc(maxnp), soln(maxnp), sumwtd, cur(maxnp)
      COMMON /jac/ ajac(3,3), resv(3)
      COMMON /pot/ phi(3), ar2
      COMMON /gstore/ gradx(maxne), grady(maxne), vsum(maxne)
