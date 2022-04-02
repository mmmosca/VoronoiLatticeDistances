#ifndef _EIGENSUPPORT_H
#define _EIGENSUPPORT_H

struct EigenVCompare {
	bool operator() (const Eigen::VectorXd& u, const Eigen::VectorXd& v) const {
		bool cond = false;
		for (int i = 0; i < u.size(); ++i) {
			if (u(i) < v(i)) {
				cond = true;
				break;
			}
			else if (u(i) > v(i)) {
				cond = false;
				break;
			}
		}
		return cond;
	}
};

#endif // !_EIGENSUPPORT_H