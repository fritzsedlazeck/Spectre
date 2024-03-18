import random
import string


class CNV_ID(object):

    def id_generator(self, size=8, chars=string.ascii_uppercase + string.digits):
        """
        Generates a random ID
        :param size: size of the ID
        :param chars: vocabulary (default only upper case chars)
        :return: random generated id with given size
        """
        return ''.join(random.choice(chars) for _ in range(size))

    @classmethod
    def n_id_generator(cls, existing_ids, n=1, size=8, chars=string.ascii_uppercase + string.digits):
        """
        Generates n amounts of unique ids with given size
        :param n: amount of new unique ids
        :param existing_ids: list of already existing ids
        :param size:  sie of the ID
        :param chars: vocabulary (default only upper case chars)
        :return: list of random generated id with given size
        """
        new_ids = []

        for i in range(0, n):
            while True:
                id = cls.id_generator(n, size, chars)
                if id not in existing_ids:
                    new_ids.append(id)
                    existing_ids.append(id)
                    break

        return new_ids
