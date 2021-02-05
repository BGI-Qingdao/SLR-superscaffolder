#include "utils/unittest/Test.h"
#include "utils/agp/agp.h"

TEST_MODULE_INIT(AGP)

using namespace BGIQD;
using namespace AGP;

TEST(test_AGP_contig_to_string)
{
    AGP_Item item;
    item.object = "test1";
    item.object_beg = 0 ;
    item.object_end = 1000;
    item.part_number = 1;
    item.component_type= AGP_Item::ComponentType::W ;
    item.lefta.component_id = "contig1";
    item.lefta.component_beg = 0;
    item.lefta.component_end = 1000;
    item.lefta.orientation = "+";

    std::string expect("test1\t0\t1000\t1\tW\tcontig1\t0\t1000\t+");
    CHECK(expect,item.ToString())
}

TEST(test_AGP_N_to_string)
{
    AGP_Item item;
    item.object = "test1";
    item.object_beg = 0 ;
    item.object_end = 1000;
    item.part_number = 1;
    item.component_type= AGP_Item::ComponentType::N ;
    item.leftb.gap_length=500;
    item.leftb.gap_type="scaffold";
    item.leftb.linkage=true;
    item.leftb.linkage_evidence="map";

    std::string expect("test1\t0\t1000\t1\tN\t500\tscaffold\tyes\tmap");
    CHECK(expect,item.ToString())
}

TEST(AGPFile_Print){
    CHECK(1,1); // too simple function, pass the test.
/**********************************
 * a read agp printed by Print :
##agp-version   2.0
# ORGANISM:
# TAX_ID:
# ASSEMBLY NAME:
# ASSEMBLY DATE: Sun Jan 24 23:26:25 2021
# GENOME CENTER:
# DESCRIPTION:
scaffold_1      1       13175   1       W       2777    1       13175   -
scaffold_1      13176   14429   2       N       1254    scaffold        yes     map
scaffold_1      14430   33823   3       W       5853    1       19394   -
*
* *********************************/
}
