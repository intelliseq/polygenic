import numpy as np
from typing import Iterable


class QuantitativeCategory(object):
    def __init__(self, category_name, from_=-np.inf, to=np.inf, scale_from = None, scale_to = None):
        super().__init__()
        self.from_ = from_
        self.to = to
        self.name = category_name
        self.scale_from = scale_from
        self.scale_to = scale_to

    def intersects_with(self, other) -> None:
        if self.from_ == other.from_ or self.to == other.to:
            return True
        elif self.from_ < other.from_ < self.to:
            return True
        elif other.from_ < self.from_ < other.to:
            return True
        return False

    def scale_cat(self, other) -> float:
        if self.scale_from is None or self.scale_to is None:
            return other
        else:
            return self.scale_from + (self.scale_to - self.scale_from) * ((other - self.from_) / (self.to - self.from_))


    def __repr__(self):
        return f'{self.__class__.__name__}({self.name}, {self.from_}, {self.to})'


class InvalidCategoriesError(RuntimeError):
    pass


def validate_categories_names(categories: Iterable[QuantitativeCategory]):
    sorted_according_to_names = sorted(categories, key=lambda x: x.name)
    for i in range(len(sorted_according_to_names) - 1):
        if categories_have_same_names(sorted_according_to_names[i], sorted_according_to_names[i + 1]):
            raise InvalidCategoriesError(
                'Categories {} and {} have the same name'.format(sorted_according_to_names[i],
                                                                 sorted_according_to_names[i + 1]))

def categories_have_same_names(category1: QuantitativeCategory, category2: QuantitativeCategory):
    return category1.name == category2.name

def validate_categories_ranges(categories: Iterable[QuantitativeCategory]):
    for i in range(len(categories) - 1):
        for j in range(i + 1, len(categories)):
            if categories[i].intersects_with(categories[j]):
                raise InvalidCategoriesError('Categories {} and {} intersect'.format(categories[i], categories[j]))