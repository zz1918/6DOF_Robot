// SSS.h

#include "SSS.cpp"

class DeltaPredicate;

class SE3Tree;

// Config is VectorXd, Box is SE3Box, BoxTree is SE3Tree, Predicate is DeltaPredicate, FeatureSet is DeltaFeature.
template<typename Config, typename Box, typename BoxTree, typename Predicate, typename FeatureSet>
class SSS;
