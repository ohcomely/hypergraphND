#include "hypergraph_ordering.hpp"

namespace hypergraph_ordering
{

    HypergraphOrdering::HypergraphOrdering(const OrderingConfig &config)
        : config_(config), stats_(), timer_()
    {
        // Validate the configuration
        if (!config.isValid())
        {
            throw std::invalid_argument("Invalid ordering configuration");
        }

        if (config_.verbose)
        {
            std::cout << "HypergraphOrdering initialized successfully." << std::endl;
            std::cout << "Configuration:" << std::endl;
            std::cout << "  Max recursion depth: " << config_.max_recursion_depth << std::endl;
            std::cout << "  Min subproblem size: " << config_.min_subproblem_size << std::endl;
            std::cout << "  Min nodes for partitioning: " << config_.min_nodes_for_partitioning << std::endl;
            std::cout << "  Imbalance tolerance: " << config_.imbalance << std::endl;
            std::cout << "  KaHyPar config: " << config_.kahypar_config_path << std::endl;
        }
    }

    kahypar_context_t *HypergraphOrdering::createKaHyParContext() const
    {
        kahypar_context_t *context = kahypar_context_new();
        kahypar_configure_context_from_file(context, config_.kahypar_config_path.c_str());

        // TODO: Can set KaHyPar seed as random
        kahypar_set_seed(context, 42);

        kahypar_supress_output(context, true);

        return context;
    }

    void HypergraphOrdering::destroyKaHyParContext(kahypar_context_t *context) const
    {
        if (context)
        {
            kahypar_context_free(context);

            if (config_.verbose)
            {
                std::cout << "KaHyPar context destroyed." << std::endl;
            }
        }
    }

    std::vector<PartitionID> HypergraphOrdering::partitionHypergraph(const HypergraphData &hg) const
    {
        Timer timer;

        // Create KaHyPar context
        kahypar_context_t *context = createKaHyParContext();

        try
        {
            if (!kahypar_validate_input(hg.num_nodes, hg.num_nets,
                                        hg.hyperedge_indices.data(), hg.hyperedges.data(),
                                        hg.net_weights.data(), hg.node_weights.data(), false))
            {
                throw std::runtime_error("Invalid hypergraph input for KaHyPar");
            }

            std::vector<PartitionID> partition(hg.num_nodes);
            kahypar_hyperedge_weight_t objective = 0;

            kahypar_partition(
                hg.num_nodes,                // num_vertices,
                hg.num_nets,                 // num_hyperedges
                config_.imbalance,           // epsilon
                config_.num_blocks,          // num_blocks
                hg.node_weights.data(),      // vertex_weights
                hg.net_weights.data(),       // hyperedge_weights
                hg.hyperedge_indices.data(), // hyperedge_indices
                hg.hyperedges.data(),        // hyperedges
                &objective,                  // objective
                context,                     // kahypar_context
                partition.data()             // partition
            );

            stats_.partitioning_time += timer.elapsed();
            stats_.hypergraph_partitions++;

            if (config_.verbose)
            {
                std::cout << "KaHyPar partitioning completed:" << std::endl;
                std::cout << "  Objective (cut): " << objective << std::endl;
                std::cout << "  Time: " << timer.elapsed() << " seconds" << std::endl;
            }

            destroyKaHyParContext(context);

            return partition;
        }
        catch (const std::exception &e)
        {
            destroyKaHyParContext(context);
            throw std::runtime_error("Error during hypergraph partitioning: " + std::string(e.what()));
        }
    }
}