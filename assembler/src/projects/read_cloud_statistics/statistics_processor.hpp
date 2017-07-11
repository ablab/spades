#pragma once

namespace read_cloud_statistics {
class Statistic {
    const string name_;

 public:
    Statistic(const string name): name_(name) {}
    virtual void Serialize(ofstream& fout) = 0;

    string GetName() const {
        return name_;
    }
};

class StatisticProcessor {
    vector<shared_ptr<Statistic>> statistics_;
    const string name_;

 public:

    StatisticProcessor(const string& name) : statistics_(), name_(name) {}

    virtual void FillStatistics() = 0;

    void SerializeStatistics(const string& base_path) {
        const string path = fs::append_path(base_path, name_);
        fs::make_dir(path);
        INFO(base_path)
        INFO(path);
        for (const auto& stat_pointer: statistics_) {
            const string new_path = fs::append_path(path, stat_pointer->GetName());
            INFO(new_path);
            ofstream stream(new_path);
            stat_pointer->Serialize(stream);
        }
    }

 protected:

    void AddStatistic(shared_ptr<Statistic> statistic) {
        statistics_.push_back(statistic);
    }

};

}