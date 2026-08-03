#include <algorithm>
#include <fstream>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>

#include <AMReX_Vector.H>
#include <gtest/gtest.h>

#include "ERF_Plotfile2DCatalog.H"
#include "ERF_Plotfile2DWaterPath.H"
#include "ERF_DataStruct.H"

using namespace plotfile2d;

namespace
{

int catalog_index (const std::string& name)
{
    const auto& catalog = diagnostic_catalog();
    for (int i = 0; i < static_cast<int>(catalog.size()); ++i) {
        if (name == catalog[i].name) {
            return i;
        }
    }
    return -1;
}

bool has_name (const amrex::Vector<std::string>& names, const std::string& target)
{
    return std::find(names.begin(), names.end(), target) != names.end();
}

// The 2D writer keeps one explicit fill block per built-in diagnostic, and the
// component index it advances is positional. The only machine-checkable link
// between the catalog and those blocks is the source text of the writer itself,
// so this test reads it. The path is derived from __FILE__ because ERF adds
// unit-test sources with absolute paths (Tests/Unit/CMakeLists.txt uses
// ${CMAKE_CURRENT_SOURCE_DIR}).
std::string plotfile2d_writer_source_path ()
{
    const std::string test_path = __FILE__;
    const std::string marker = "Tests/Unit/IO/";
    const auto pos = test_path.rfind(marker);
    if (pos == std::string::npos) {
        return {};
    }
    return test_path.substr(0, pos) + "Source/IO/ERF_Plotfile2D.cpp";
}

// Returns the built-in diagnostic names in the order their fill blocks appear
// in the 2D writer. Land-surface fields are filled by a generic name-driven loop
// and deliberately have no literal here.
std::vector<std::string> writer_fill_order (const std::string& source_text)
{
    const std::string needle = "containerHasElement(plot_var_names, \"";
    std::vector<std::string> names;

    std::size_t pos = source_text.find(needle);
    while (pos != std::string::npos) {
        const std::size_t name_begin = pos + needle.size();
        const std::size_t name_end = source_text.find('"', name_begin);
        if (name_end == std::string::npos) {
            break;
        }
        names.push_back(source_text.substr(name_begin, name_end - name_begin));
        pos = source_text.find(needle, name_end);
    }

    return names;
}

} // namespace

// Motivation: urb_frac is the tile blend weight the urban scheme will key off,
// so it must be a first-class geometry field with stable metadata rather than a
// scheme-conditional land-surface diagnostic.
TEST(UrbanFraction, CatalogRegistersUrbFracAsGeometry)
{
    const auto* descriptor = plotfile2d::find_diagnostic("urb_frac");
    ASSERT_NE(descriptor, nullptr);

    EXPECT_EQ(descriptor->id, plotfile2d::DiagnosticID::UrbFrac);
    EXPECT_EQ(descriptor->category, plotfile2d::DiagnosticCategory::Geometry);
    EXPECT_EQ(descriptor->missing_policy, plotfile2d::MissingPolicy::AlwaysAvailable);
    EXPECT_STREQ(descriptor->units, "1");
    EXPECT_FALSE(std::string(descriptor->long_name).empty());
}

// Motivation: The catalog order is the canonical plotfile component order, and
// the writer fills urb_frac immediately after landmask. Pinning the adjacency
// here makes an accidental catalog reorder fail in the unit suite instead of
// silently shifting every later component in written plotfiles.
TEST(UrbanFraction, UrbFracFollowsLandmaskInCatalog)
{
    const int landmask_index = catalog_index("landmask");
    const int urb_frac_index = catalog_index("urb_frac");
    const int mapfac_index = catalog_index("mapfac");

    ASSERT_GE(landmask_index, 0);
    ASSERT_GE(urb_frac_index, 0);
    ASSERT_GE(mapfac_index, 0);

    EXPECT_EQ(urb_frac_index, landmask_index + 1);
    EXPECT_EQ(mapfac_index, urb_frac_index + 1);
}

// Motivation: urb_frac_lev is allocated for every level regardless of the land
// model, so the field must be selectable in every configuration. A container
// that is always allocated but conditionally advertised would make the
// requested-vs-written component counts diverge.
TEST(UrbanFraction, UrbFracIsAvailableInEveryConfiguration)
{
    SolverChoice dry;
    dry.moisture_type = MoistureType::None;
    EXPECT_TRUE(has_name(plotfile2d::available_diagnostic_names(dry), "urb_frac"));
    EXPECT_TRUE(has_name(plotfile2d::available_diagnostic_names(dry, false), "urb_frac"));
    EXPECT_TRUE(has_name(plotfile2d::available_diagnostic_names(dry, true), "urb_frac"));

    SolverChoice moist;
    moist.moisture_type = MoistureType::MoistNoCondensation;
    const amrex::Vector<std::string> active_lsm{"t_sfc", "smois_1"};
    EXPECT_TRUE(has_name(plotfile2d::available_diagnostic_names(moist, true, active_lsm),
                         "urb_frac"));
}

// Motivation: Catalog order and writer fill order are coupled only by
// convention; when they disagree the writer either aborts on the component
// count or, worse, labels one field with another field's name. This is the
// regression assertion for that coupling.
TEST(UrbanFraction, WriterFillOrderMatchesCatalogOrder)
{
    const std::string source_path = plotfile2d_writer_source_path();
    ASSERT_FALSE(source_path.empty())
        << "Could not derive the 2D writer path from __FILE__ = " << __FILE__;

    std::ifstream source_file(source_path);
    ASSERT_TRUE(source_file.good())
        << "Could not open the 2D writer source at " << source_path
        << ". This test reads it to compare fill order against the catalog.";

    std::stringstream buffer;
    buffer << source_file.rdbuf();
    const std::vector<std::string> fill_order = writer_fill_order(buffer.str());
    ASSERT_FALSE(fill_order.empty())
        << "Found no fill blocks in " << source_path
        << "; the containerHasElement(plot_var_names, \"...\") idiom may have changed.";

    // Every fill block must name a catalog entry, and the blocks must appear in
    // strictly increasing catalog order.
    int previous_index = -1;
    std::string previous_name = "<start of writer>";
    for (const auto& name : fill_order) {
        const int index = catalog_index(name);
        EXPECT_GE(index, 0) << "Writer fills '" << name << "', which is not in the catalog";
        if (index < 0) {
            continue;
        }
        EXPECT_GT(index, previous_index)
            << "Writer fills '" << name << "' after '" << previous_name
            << "', but the catalog orders them the other way";
        previous_index = index;
        previous_name = name;
    }

    // No diagnostic may be filled twice: each fill block advances mf_comp once.
    std::vector<std::string> sorted_fill_order = fill_order;
    std::sort(sorted_fill_order.begin(), sorted_fill_order.end());
    EXPECT_EQ(std::adjacent_find(sorted_fill_order.begin(), sorted_fill_order.end()),
              sorted_fill_order.end())
        << "A 2D diagnostic has more than one fill block";

    // Every built-in that is not filled by the generic land-surface loop or by
    // the dynamic sampled-level path needs an explicit fill block. This is what
    // catches a catalog entry added without a corresponding writer block.
    for (const auto& descriptor : diagnostic_catalog()) {
        if (descriptor.category == plotfile2d::DiagnosticCategory::LandSurface ||
            descriptor.category == plotfile2d::DiagnosticCategory::SampledLevel) {
            continue;
        }
        EXPECT_NE(std::find(fill_order.begin(), fill_order.end(), descriptor.name),
                  fill_order.end())
            << "Catalog entry '" << descriptor.name
            << "' has no fill block in " << source_path;
    }

    // And specifically: urb_frac is filled, and it is filled between landmask
    // and mapfac.
    const auto landmask_it = std::find(fill_order.begin(), fill_order.end(), "landmask");
    const auto urb_frac_it = std::find(fill_order.begin(), fill_order.end(), "urb_frac");
    const auto mapfac_it = std::find(fill_order.begin(), fill_order.end(), "mapfac");
    ASSERT_NE(landmask_it, fill_order.end());
    ASSERT_NE(urb_frac_it, fill_order.end());
    ASSERT_NE(mapfac_it, fill_order.end());
    EXPECT_EQ(std::distance(landmask_it, urb_frac_it), 1);
    EXPECT_EQ(std::distance(urb_frac_it, mapfac_it), 1);
}
