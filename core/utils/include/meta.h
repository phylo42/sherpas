#ifndef RAPPAS_CORE_META_H
#define RAPPAS_CORE_META_H

namespace rappas
{
    template <bool flag, class IsTrue, class IsFalse>
    struct choose;

    template <class IsTrue, class IsFalse>
    struct choose<false, IsTrue, IsFalse> {
        typedef IsFalse type;
    };

    template <class IsTrue, class IsFalse>
    struct choose<true, IsTrue, IsFalse> {
        typedef IsTrue type;
    };
}

#endif //RAPPAS_CORE_META_H
